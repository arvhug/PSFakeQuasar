import numpy as np
from scipy.stats import lognorm
from astropy.table import Table
from astropy.stats import SigmaClip
from photutils.psf import PSFPhotometry

def psfInjection(fake_flux, psf, psf_model, ratio, cutout, x0, y0, ZP=None):
    """
    Inject a PSF into a cutout image at specified coordinates.
    
    Parameters:
    -----------
    fake_flux : float
        Desired flux of the injected PSF (in microJanskys if ZP is provided, otherwise in ADU)
    psf : ndarray
        The PSF image data
    psf_model : FittableImageModel
        The PSF model to evaluate
    ratio : float
        Oversampling ratio of the PSF
    cutout : Cutout2D or ndarray
        The image cutout where PSF will be injected
    x0, y0 : float
        Central coordinates where PSF should be injected
    ZP : float, optional
        Zero point magnitude for flux conversion. If provided, fake_flux is assumed to be in μJy.
    
    Returns:
    --------
    ndarray
        The cutout image with injected PSF
    """
    # Initialize fake PSF image with same shape as cutout
    if hasattr(cutout, 'data'):  # Handle both Cutout2D and raw array inputs
        cutout_data = cutout.data
    else:
        cutout_data = cutout
    fake_ps = np.zeros_like(cutout_data)
    
    # Calculate PSF size and position
    size_psf = psf.shape[0] // int(ratio)
    x_add = x0 - size_psf//2
    y_add = y0 - size_psf//2
    
    # Check if PSF will fit within image bounds
    if not (x_add >= 0 and y_add >= 0 and 
            x_add + size_psf <= fake_ps.shape[1] and 
            y_add + size_psf <= fake_ps.shape[0]):
        raise ValueError("PSF would be placed outside the image bounds")
    
    # Create coordinate grid for PSF evaluation
    xx, yy = np.meshgrid(np.arange(size_psf), np.arange(size_psf))
    xc_psf, yc_psf = size_psf//2, size_psf//2
    
    # Adjust center for even-sized PSF
    if size_psf % 2 == 0:
        xc_psf = xc_psf - 0.5
        yc_psf = yc_psf - 0.5
    
    #Generate random offset (astrometric uncertainty)
    shape, loc, scale = 0.4165984893759021, -0.2698935690678829, 0.8961365668573444
    offset = np.random.uniform(-0.5, 0.5, 2)
    distances = lognorm(s=shape, loc=loc, scale=scale).rvs(1)[0]
    angles = np.random.uniform(0, 2 * np.pi, 1)[0]
    dx = distances * np.cos(angles)
    dy = distances * np.sin(angles)
    xc_psf += dx
    yc_psf += dy
    
    # Evaluate PSF model
    _psf = psf_model.evaluate(xx, yy, flux=1, x_0=xc_psf+offset[0], y_0=yc_psf+offset[1], use_oversampling=True)
    #_psf = psf_model.evaluate(xx, yy, flux=1, x_0=xc_psf, y_0=yc_psf, use_oversampling=True)
    
    # Place PSF in the image
    fake_ps[y_add:y_add+size_psf, x_add:x_add+size_psf] = _psf
    
    # Convert flux if ZP is provided
    if ZP is not None:
        # Convert μJy to ADU using zero point
        fake_flux =fake_flux * 10**(0.4 * (ZP - 23.9))  
    
    # Scale PSF to desired flux
    fake_ps = fake_ps * fake_flux / np.sum(fake_ps)
    
    # Add PSF to cutout
    _img = cutout_data + fake_ps
    
    return _img

def generate_random_positions_and_psf_flux(psf_phot, source_x, source_y, 
                                         inner_radius_arcseco, outer_radius_arcseco, 
                                         num_apertures, data, 
                                         pixel_scale):
    """
    Generate random positions around the source and perform PSF photometry
    to estimate background fluctuations.
    
    Parameters:
    -----------
    psf_model : Astropy PSF model
        The PSF model to use for photometry
    source_x, source_y : float
        Source coordinates in pixel coordinates
    inner_radius_arcsec, outer_radius_arcsec : float
        Annulus for random position generation (in arcsec)
    num_apertures : int
        Number of random positions to generate
    cutout_data : ndarray
        Image data
    back_mean : float
        Background mean value to subtract
    pixel_scale : float
        Pixel scale in arcsec/pixel
    aperture_radius : float
        Aperture radius for PSF photometry
        
    Returns:
    --------
    tuple: (list of flux measurements, list of flux error measurements)
    """
    flux_measurements = []
    flux_error_measurements = []
    
    # Convert radii from arcsec to pixels
    inner_radius = inner_radius_arcseco / pixel_scale
    outer_radius = outer_radius_arcseco / pixel_scale
    
    # Set up PSF photometry (same as in your main function)
    # psf_model.x_0.fixed = False
    # psf_model.y_0.fixed = False
    # psf_phot = PSFPhotometry(psf_model, (21, 21), aperture_radius=aperture_radius)
    
   
    
    for i in range(num_apertures):
        # Generate random position in annulus
        r = np.sqrt(np.random.uniform(inner_radius**2, outer_radius**2))
        theta = np.random.uniform(0, 2 * np.pi)
        
        x_rand = source_x + r * np.cos(theta)
        y_rand = source_y + r * np.sin(theta)
        
        # Check if position is within image bounds (with some buffer for PSF)
        buffer = 10  # pixels, adjust based on PSF size
        if (x_rand < buffer or x_rand > data.shape[1] - buffer or
            y_rand < buffer or y_rand > data.shape[0] - buffer):
            continue
            
        # Perform PSF photometry at random position
        init_params = Table()
        init_params['x_0'] = [x_rand]
        init_params['y_0'] = [y_rand]
        
        try:
            phot = psf_phot(data=data, init_params=init_params)
            flux_psf = phot['flux_fit'][0]
            flux_psf_err = phot['flux_err'][0]
            
            flux_measurements.append(flux_psf)
            flux_error_measurements.append(flux_psf_err)
        except:
            continue
    
    return flux_measurements, flux_error_measurements

def calculate_psf_flux_vars(flux_measurements, zero_point):
    """
    Calculate statistics from PSF flux measurements.
    
    Parameters:
    -----------
    flux_measurements : list
        List of flux measurements from background positions
    zero_point : float
        Zero point magnitude
        
    Returns:
    --------
    tuple: (clipped_std_flux, median_flux_err) in microJanskys
    """
    if not flux_measurements:
        return np.nan, np.nan
    
    # Sigma clip to remove outliers
    sigma_clip_i = SigmaClip(sigma=3, maxiters=10)
    clipped_flux = sigma_clip_i(flux_measurements)
    
    # Calculate statistics
    std_flux = np.std(clipped_flux)
    std_flux_uJy = std_flux * np.power(10, -0.4 * (zero_point - 23.9))
    
    return std_flux_uJy


def perform_psf_photometry(psf_model, xc, yc, cutout, back_mean, 
                         pixel_scale, seg=None, ZP=23.9, aperture_radius=4,
                         bg_annulus=(2, 5), num_bg_samples=30):
    """
    Updated version with PSF-based background error estimation.
    bg_annulus: tuple (inner, outer radius in arcsec)
    num_bg_samples: number of background positions to sample
    """
    # Main PSF photometry (unchanged from your original)
    psf_model.x_0.fixed = False
    psf_model.y_0.fixed = False
    fitshape = (21, 21)
    psf_phot = PSFPhotometry(psf_model, fit_shape=fitshape, aperture_radius=aperture_radius)
    init_params = Table()
    init_params['x_0'] = [xc]
    init_params['y_0'] = [yc]

    data_to_fit = cutout.data - back_mean

    try:
        phot = psf_phot(data=data_to_fit, mask=seg, init_params=init_params)
        flux_psf = phot['flux_fit'][0]
    except:
        flux_psf = np.nan
        flux_psf_err = np.nan
        flux_psf_uJy = np.nan
        flux_psf_err_uJy = np.nan
        return flux_psf_uJy, flux_psf_err_uJy, phot
    
    # PSF-based background error estimation
    bg_fluxes, bg_errors = generate_random_positions_and_psf_flux(
        psf_phot, xc, yc, 
        bg_annulus[0], bg_annulus[1],
        num_bg_samples, data_to_fit, 
        pixel_scale
    )

    #print(bg_fluxes)
    
    if bg_fluxes:
        flux_err = calculate_psf_flux_vars(bg_fluxes, ZP)
    else:
        flux_err = np.nan
    
    # Convert to microJanskys
    flux_psf_uJy = flux_psf* 10**(-0.4*(ZP-23.9))
    flux_psf_err_uJy = flux_err  # Already in uJy from calculate_psf_flux_vars

    mag_psf = float(-2.5 * np.log10(flux_psf) + ZP)
    
    return flux_psf_uJy, flux_psf_err_uJy, phot, mag_psf
