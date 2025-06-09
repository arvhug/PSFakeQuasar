# PSFakeQuasar :comet: :syringe: :stars: :dart: 

> **"Not all quasars are real... but yours can be!"**  


Tools for PSF modelling, injection, and photometry.
Tailored to the Euclid survey.

## Local Installation
```bash
git clone https://github.com/arvhug/PSFakeQuasar.git
cd PSFakeQuasar
pip install -e .
python -c "from PSFakeQuasar import inject_quasar; print('Success!')"
```


## Usage

```python
from PSFakeQuasar import inject_quasar, load_psf

# Load a Euclid PSF
psf = load_psf("euclid_psf.fits", ra=53.2, dec=-27.8)  

# Inject a 100 μJy quasar at (x=50, y=60)
image = inject_quasar(
    fake_flux=100, 
    psf=psf, 
    x0=50, y0=60, 
    ZP=24.5
)
```
## License
This project is licensed under the **GNU GPLv3**. See [LICENSE](LICENSE) for details.