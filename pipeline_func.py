from astropy.coordinates import Angle
from astropy import units as u
from astropy.io import fits
import numpy as np
from collections import defaultdict

def tree():
    return defaultdict(tree)

# Process the FITS files for all objects
def process_obj_astrometry(info):

    # Location of FITS in memory
    rootdir = info[0] 
    # Object name
    obj_name = info[1]
    # Pool queue
    q = info[2]
    
    # Load object data
    obj_info = tree()

    # Open the FITS file
    obj_fits = fits.open(rootdir + "Spectra/" + obj_name + ".fits")
    
    # Obtain RA, Dec
    primary_hdr = obj_fits[0].header
    ra = Angle(primary_hdr['plug_ra']*u.deg).hms
    dec = Angle(primary_hdr['plug_dec']*u.deg).dms
    
    obj_info['ra'] = ra
    obj_info['dec'] = dec
    
    # Obtain redshift
    spec_obj = obj_fits[2].data
    spec_hdr = obj_fits[2].header

    n_fields = spec_hdr['tfields']

    redshift_idx = 0
 
    for i in range(1, n_fields):
        if(spec_hdr['ttype' + str(i)] == 'Z'):
            redshift_idx = i
            break
    
    redshift = spec_obj[0][redshift_idx-1]
    obj_info['redshift'] = redshift
 
    #print('Completed rendering object {} \n'.format(str(obj_name)))

    q.put((obj_name, obj_info))

def emission_pipeline(info):

    # Unpack info
    rootdir = info[0]
    obj_name = info[1]
    q = info[2]
    line_ids = info[3]

    # Load FITS file
    obj_fits = fits.open(rootdir + "Spectra/" + obj_name + ".fits")
    emission_info = obj_fits[3].data

    # Locate where the line is
    n_lines = len(emission_info)
    
    line_info = tree()
    
    for i in range(0, n_lines):
        if(emission_info[i][3] in line_ids):
            # Extract flux and EW
            line_info[emission_info[i][3]] = [emission_info[i][9], emission_info[i][11]]
    
    #print('Completed rendering object {} \n'.format(str(obj_name)))

    q.put((obj_name, line_info))

def emission_custom(info):

    # Unpack info
    rootdir = info[0]
    obj_name = info[1]
    q = info[2]
    line_windows = info[3]

    # Open the FITS file
    obj_fits = fits.open(rootdir + "Spectra/" + obj_name + ".fits")

    #Obtain spectroscopic information
    spec = list(obj_fits[1].data)
    spec = [list(spec[i]) for i in range(len(spec))]
    spec = np.array(spec)
    
    wav = 10**spec[:,1]
    flux = spec[:,0]

    # Obtain redshift
    spec_obj = obj_fits[2].data
    spec_hdr = obj_fits[2].header

    n_fields = spec_hdr['tfields']

    redshift_idx = 0
 
    for i in range(1, n_fields):
        if(spec_hdr['ttype' + str(i)] == 'Z'):
            redshift_idx = i
            break
    
    redshift = spec_obj[0][redshift_idx-1]

    # Rest frame wavelength
    wav_shifted = wav/(1+redshift)

    # Wavelength pixel width
    wav_del = np.average(wav_shifted[1:]-wav_shifted[:-1])

    # Get list of lines
    line_ids = list(line_windows.keys())

    line_info = tree()

    # Extract line windows
    for l in line_ids:

        # Central wavelength
        cent = line_windows[l][0]

        # Width of left wing
        left_width = line_windows[l][1]
        # Left-wing to left continuum
        left_space = line_windows[l][2]
        # Left continuum width
        left_cont = line_windows[l][3]

        # Width of right wing
        right_width = line_windows[l][4]
        # Right-wing to right continuum
        right_space = line_windows[l][5]
        # Right continuum width
        right_cont = line_windows[l][6]

        # Left most wavelength
        min_wav = cent - left_width - left_space - left_cont
        # Right most wavelength
        max_wav = cent + right_width + right_space + right_cont

        if(wav_shifted[0]>min_wav or wav_shifted[-1]<max_wav):
            line_info[l] = [0., 0., 0., 0.]

        else:
            # Get continuum indices
            left_cont_idx = ((wav_shifted>=min_wav) 
                                & (wav_shifted<=(min_wav+left_cont)))

            right_cont_idx = ((wav_shifted>=(max_wav-right_cont)) 
                                & (wav_shifted<=max_wav))

            # Get emission wing indices
            i = ((wav_shifted>=(cent-left_width)) & (wav_shifted<=(cent+right_width)))

            # Central pixel index
            central_idx = (np.abs(wav_shifted - cent)).argmin()

            # Calculate line intensity
            wav_line = wav_shifted[i]
            flux_line = flux[i]

            line_intensity = np.trapz(flux_line, x=wav_line, dx=wav_del)

            # Fit straight line to continuum

            left_wav_cont = wav_shifted[left_cont_idx]
            right_wav_cont = wav_shifted[right_cont_idx]

            mid_left_wav = left_wav_cont[len(left_wav_cont)//2]
            mid_right_wav = right_wav_cont[len(right_wav_cont)//2]

            left_flux_cont = flux[left_cont_idx]
            right_flux_cont = flux[right_cont_idx]

            med_left_flux = np.median(left_flux_cont)
            med_right_flux = np.median(right_flux_cont)

            slope = (med_right_flux-med_left_flux)/(mid_right_wav-mid_left_wav)
            intercept = med_left_flux - slope*mid_left_wav

            # Compute continuum level at central wavelength, continuum intensity, alongside equivalent width
            cont_flux = slope*wav_shifted[central_idx] + intercept
            cont_intensity = np.trapz(slope*wav_line + intercept, x=wav_line, dx=wav_del)

            # Line equivalent width
            ew = (line_intensity-cont_intensity)/cont_flux

            # Continuum noise and SNR
            wav_cont = np.concatenate((left_wav_cont, right_wav_cont))
            flux_cont = np.concatenate((left_flux_cont, right_flux_cont))

            noise = np.sqrt(np.average((flux_cont - (slope*wav_cont + intercept))**2))
            snr_c = cont_flux/noise

            # Line SNR
            snr_l = (line_intensity-cont_intensity)/(noise*wav_del*np.sqrt(len(wav_line)))

            line_info[l] = [line_intensity-cont_intensity, ew, snr_l, snr_c]

        #print('Completed rendering object {} \n'.format(str(obj_name)))

    q.put((obj_name, line_info))