import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
from astroquery.vizier import Vizier
import subprocess
import os
import sys

def test_dependency(dep, alternate_name=None):
    """
    Test external dependency by trying to run it as a subprocess
    """
    try:
        subprocess.check_output(dep, stderr=subprocess.PIPE, shell=True)
        print("%s is installed properly as %s. OK" % (dep, dep))
        return 1
    except subprocess.CalledProcessError:
        try:
            subprocess.check_output(alternate_name, stderr=subprocess.PIPE, shell=True)
            print("%s is installed properly as %s. OK" % (dep, alternate_name))
            return 1
        except subprocess.CalledProcessError:
            print("===%s/%s IS NOT YET INSTALLED PROPERLY===" % (dep, alternate_name))
            return 0
    
dependencies = [('source-extractor', 'sex'), ('psfex', 'PSFEx')]
i = 0
for dep_name1, dep_name2 in dependencies:
    i += test_dependency(dep_name1, dep_name2)
print("%i out of %i external dependencies installed properly.\n" % (i, len(dependencies)))
if i != len(dependencies):
    print("Please correctly install these programs before continuing by following the instructions in README.md.")
else:
    print("You are ready to continue.") 

os.chdir('/home/aditi/ksp/Task_2/files_TASK_2/')
imageName = sys.argv[1]

with fits.open(imageName) as HDUList:
    header = HDUList[0].header
    image = HDUList[0].data

zscale = ZScaleInterval().get_limits(image)

plt.figure(figsize=(10,10))
plt.imshow(image, cmap='gray', origin='lower', vmin=zscale[0], vmax=zscale[1])
plt.colorbar()
plt.show()

w = WCS(header)
(raImage, decImage) = w.all_pix2world(image.shape[0]/2, image.shape[1]/2, 1)
boxsize = 30 # arcminutes
maxmag = 18

catNum = 'II/349'
print(f'\nQuerying Vizier {catNum} around RA {raImage:.4f}, Dec {decImage:.4f} with a radius of {boxsize} arcmin')

try:
    # You can set the filters for the individual columns (magnitude range, number of detections) inside the Vizier query
    v = Vizier(columns=['*'], column_filters={"gmag":f"<{maxmag}", "Nd":">6", "e_gmag":f"<{1.086/3}"}, row_limit=-1)
    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), radius=str(boxsize)+'m', catalog=catNum, cache=False)
    print(Q[0])
except:
    print('I cannnot reach the Vizier database. Is the internet working?')

Q[0].meta['desc'] = Q[0].meta.pop('description')
Q[0].write('/home/aditi/ksp/Task_2/files_TASK_2/ps1Catalog.fits', format='fits', overwrite=True)

ps1_imCoords = w.all_world2pix(Q[0]['RAJ2000'], Q[0]['DEJ2000'], 1)
good_cat_stars = Q[0][np.where((ps1_imCoords[0] > 500) & (ps1_imCoords[0] < 3500) & (ps1_imCoords[1] > 500) & (ps1_imCoords[1] < 3500))]
ps1_imCoords = w.all_world2pix(good_cat_stars['RAJ2000'],good_cat_stars['DEJ2000'], 1)

fig = plt.figure(figsize=(10,10))
ax = fig.gca()
plt.imshow(image, cmap='gray', vmin=zscale[0], vmax=zscale[1])
circles = [plt.Circle((ps1_imCoords[0][i], ps1_imCoords[1][i]), radius = 10, edgecolor='C1', facecolor='None') for i in range(len(ps1_imCoords[0]))]
for c in circles:
    ax.add_artist(c)
    
plt.show()

configFile = "/home/aditi/ksp/Task_2/files_TASK_2/photomCat.sex"
catalogName = imageName + '.cat'
paramName = "/home/aditi/ksp/Task_2/files_TASK_2/photomCat.param"
try:
    command = f'source-extractor -c {configFile} {imageName} -CATALOG_NAME {catalogName} -PARAMETERS_NAME {paramName}'
    print(f'Executing command: {command}')
    rval = subprocess.run(command.split(), check=True)
except subprocess.CalledProcessError as err:
    print(f'Could not run sextractor with exit error {err}')

with fits.open(catalogName) as HDU:
    print(HDU.info())
    sourceTable = Table(HDU[2].data)

print(sourceTable.colnames)
print(sourceTable)

cleanSources = sourceTable[(sourceTable['FLAGS']==0) & (sourceTable['FWHM_WORLD'] < 2) & (sourceTable['XWIN_IMAGE']<3500) & (sourceTable['XWIN_IMAGE']>500) & (sourceTable['YWIN_IMAGE']<3500) & (sourceTable['YWIN_IMAGE']>500)]




fig = plt.figure(figsize=(10,10))
ax = fig.gca()
plt.imshow(image, cmap='gray', vmin=zscale[0], vmax=zscale[1])
circles = [plt.Circle((source['XWIN_IMAGE'], source['YWIN_IMAGE']), radius = 10, edgecolor='C1', facecolor='None') for source in cleanSources]
for c in circles:
    ax.add_artist(c)
plt.show()


psfConfigFile = 'psfex_conf.psfex'

try:
    command = f'psfex -c {psfConfigFile} {catalogName}'
    print(f'Executing command: {command}')
    rval = subprocess.run(command.split(), check=True)
except subprocess.CalledProcessError as err:
    print(f'Could not run psfex with exit error {err}')



psfModelHDU = fits.open('moffat_' + imageName + '.fits')[0]
psfModelData = psfModelHDU.data
mean, median, std = sigma_clipped_stats(psfModelData)

plt.figure(figsize=(6,6))
plt.imshow(psfModelData, vmin=0, vmax=median+20*std, origin='lower')
plt.show()

psfName = imageName + '.psf'
psfcatalogName = imageName+'.psf.cat'
psfparamName = 'photomPSF.param' # This is a new set of parameters to be obtained from SExtractor, including PSF-fit magnitudes
try:
    # We are supplying SExtactor with the PSF model with the PSF_NAME option
    command = f'source-extractor -c {configFile} {imageName} -CATALOG_NAME {psfcatalogName} -PSF_NAME {psfName} -PARAMETERS_NAME {psfparamName}'
    print(f"Executing command: {command}")
    rval = subprocess.run(command.split(), check=True)
except subprocess.CalledProcessError as err:
    print(f'Could not run sextractor with exit error {err}')


with fits.open(psfcatalogName) as HDU:
    psfsourceTable = Table(HDU[2].data)

print(psfsourceTable.colnames)
print(psfsourceTable)

cleanPSFSources = psfsourceTable[(psfsourceTable['FLAGS']==0) & (psfsourceTable['FLAGS_MODEL']==0)  & (psfsourceTable['FWHM_WORLD'] < 2) & (psfsourceTable['XMODEL_IMAGE']<3500) & (psfsourceTable['XMODEL_IMAGE']>500) &(psfsourceTable['YMODEL_IMAGE']<3500) &(psfsourceTable['YMODEL_IMAGE']>500)]

psfsourceCatCoords = SkyCoord(ra=cleanPSFSources['ALPHAWIN_J2000'], dec=cleanPSFSources['DELTAWIN_J2000'], frame='icrs', unit='degree')
ps1CatCoords = SkyCoord(ra=good_cat_stars['RAJ2000'], dec=good_cat_stars['DEJ2000'], frame='icrs', unit='degree')
photoDistThresh = 0.6
idx_psfimage, idx_psfps1, d2d, d3d = ps1CatCoords.search_around_sky(psfsourceCatCoords, photoDistThresh*u.arcsec)

print(f'Found {len(idx_psfimage)} good cross-matches')

plt.figure(figsize=(8,8))
plt.scatter(good_cat_stars['gmag'][idx_psfps1], cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage], color='C0', s=10)
plt.xlabel('PS1 magnitude', fontsize=15)
plt.ylabel('Instrumental PSF-fit magnitude', fontsize=15)
plt.show()

psfoffsets = ma.array(good_cat_stars['gmag'][idx_psfps1] - cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage])
zero_psfmean, zero_psfmed, zero_psfstd = sigma_clipped_stats(psfoffsets)
print('PSF Mean ZP: %.2f\nPSF Median ZP: %.2f\nPSF STD ZP: %.2f'%(zero_psfmean, zero_psfmed, zero_psfstd))

ra = 210.910674637
dec = 54.3116510708

sn2023ixf_coords = SkyCoord(ra=[ra], dec=[dec], frame='icrs', unit='degree')
idx_sn2023ixf, idx_cleanpsf_sn2023ixf, d2d, d3d = psfsourceCatCoords.search_around_sky(sn2023ixf_coords, photoDistThresh*u.arcsec)
print(f'Found the source at index {idx_cleanpsf_sn2023ixf[0]}')

sn2023ixf_psfinstmag = cleanPSFSources[idx_cleanpsf_sn2023ixf]['MAG_POINTSOURCE'][0]
sn2023ixf_psfinstmagerr = cleanPSFSources[idx_cleanpsf_sn2023ixf]['MAGERR_POINTSOURCE'][0]

sn2023ixf_psfmag = zero_psfmed + sn2023ixf_psfinstmag
sn2023ixf_psfmagerr = np.sqrt(sn2023ixf_psfinstmagerr**2 + zero_psfstd**2)

print(f'PSF-fit magnitude of SN2023ixf is {sn2023ixf_psfmag} +/- {sn2023ixf_psfmagerr}')
os.chdir("/home/aditi/ksp/Task_2")
with open('data.txt', 'a') as file:
    lis=f'{imageName},{sn2023ixf_psfmag},{sn2023ixf_psfmagerr}\n'
    file.write(lis)
    print("done")