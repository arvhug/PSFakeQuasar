import numpy as np
from astropy.io import fits
from photutils.psf import FittableImageModel

# def load_psf(file_path_psf, ra_target, dec_target):
#     """Load PSF from file and extract stamp at target coordinates."""
#     from MER_PsfMosaic.MerCatalogPsf import MerCatalogPsf
#     from MER_PsfMosaic.EuclidWcs import EuclidWcs
    
#     psf = MerCatalogPsf.from_file(file_path_psf)
#     psf.set_wcs(EuclidWcs(psf.get_header()))
    
#     psf_stamp = psf.get_closest_stamp_at_radec(np.array([ra_target, dec_target]))
#     return psf_stamp.get_data()

def normPSF(psf, pixel_scale, band):
    """Normalize PSF for given pixel scale and band."""
    psf_nisp_pixscale = 0.1  # NISP pixel scale
    psf_vis_pixscale = 0.1   # VIS pixel scale
    
    psf_pixscale = psf_vis_pixscale if band == 'VIS' else psf_nisp_pixscale
    ratio = pixel_scale / psf_pixscale
    psf_norm = psf * (int(ratio)**2) / np.sum(psf)
    return FittableImageModel(psf_norm, normalize=False, oversampling=int(ratio)), ratio

def get_psf_stamp(file_path, stamp_idx=None, normalize=True):
    """
    Extract a PSF stamp from a FITS file, either random or specified.
    
    Parameters
    ----------
    file_path : str
        Path to the FITS file containing PSF stamps
    stamp_idx : int, optional
        Index of specific stamp to extract. If None, picks randomly.
    normalize : bool, optional
        Whether to normalize the stamp (sum to 1)
        
    Returns
    -------
    tuple
        (stamp_data, stamp_header, metadata) where metadata contains:
        - stamp_position: (x_start, y_start)
        - stamp_size: int
        - total_stamps: int
    """
    with fits.open(file_path) as hdul:
        psf_data = hdul[1].data
        header = hdul[1].header
    
    stamp_size = header['STMPSIZE']
    n_stamps_x = psf_data.shape[1] // stamp_size
    n_stamps_y = psf_data.shape[0] // stamp_size
    total_stamps = n_stamps_x * n_stamps_y
    
    # Handle stamp selection
    if stamp_idx is None:
        stamp_idx = np.random.randint(0, total_stamps)
    elif stamp_idx >= total_stamps:
        raise ValueError(f"Stamp index {stamp_idx} exceeds total stamps ({total_stamps})")
    
    # Calculate position
    x_start = (stamp_idx % n_stamps_x) * stamp_size
    y_start = (stamp_idx // n_stamps_x) * stamp_size
    
    # Extract stamp
    stamp = psf_data[
        y_start : y_start + stamp_size,
        x_start : x_start + stamp_size
    ]
    
    if normalize:
        stamp = stamp / stamp.sum()
    
    # Prepare metadata
    metadata = {
        'stamp_position': (x_start, y_start),
        'stamp_size': stamp_size,
        'total_stamps': total_stamps,
        'stamp_idx': stamp_idx
    }
    
    return stamp, header, metadata
