# PSFakeQuasar :comet: :syringe: :stars: :dart: 

> **"Source Injector and assorted PSF tools"**  
Tailored to the Euclid survey will generalise in the future.

## Local Installation
```bash
git clone https://github.com/arvhug/PSFakeQuasar.git
cd PSFakeQuasar
pip install -e .
python -c "from PSFakeQuasar import psfInjection; print('Success!')"
```


## Usage

Data required: Science image, PSF, Flux of Source

```python
from PSFakeQuasar import psfInjection, get_psf_stamp,normPSF


## Your science image; then cutot at a target position with a particular size
cutout = Cutout2D(sci_data, coord_target, size, wcs=wcs, mode='partial')
## entre of the image to inject
xc, yc = cutout.data.shape[1] // 2, cutout.data.shape[0] // 2
######

file_path_psf = 'your_Euclid_psf.fits
psf_data, _, _ = get_psf_stamp(file_path_psf)
psf_norm,ratio = normPSF(psf_data,'NISP',0.1)

image= psfInjection(1.9063500394148085, psf_data, psf_norm, ratio, cutout, xc, yc, ZP=ZP)

## Check the source injection
plt.figure(figsize=(10, 8))
plt.imshow(image, origin='lower', vmin=0, vmax=15, cmap='viridis')
plt.colorbar(label='Flux')


```
## License
This project is licensed under the **GNU GPLv3**. See [LICENSE](LICENSE) for details.