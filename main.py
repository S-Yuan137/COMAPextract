import COMANCextract
import COMANCextract.utils as utils
from COMANCextract.class_map import AstroMap
from COMANCextract.ap_photo import apPhotometry
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from matplotlib.backends.backend_pdf import PdfPages
import sys
from tqdm import tqdm
import h5py

utils.print_header(version=COMANCextract.__version__)
params = utils.parse_file(sys.argv[1])

for SourceID in list(params['sources'].keys()):
    new_wcs = WCS(naxis=2)
    new_wcs.wcs.crpix = np.asarray(params['wcs']['crpix']).astype(int)
    new_wcs.wcs.cdelt = np.asarray(params['wcs']['cdelt'])/60
    new_wcs.wcs.crval = np.asarray(params['sources'][SourceID]['coords'])
    new_wcs.wcs.ctype = params['wcs']['ctype']
    new_shape = (new_wcs.wcs.crpix*2 - np.ones(2)).astype(int)

    pdf = PdfPages('{}_PhotometryReport_{}.pdf'.format(
              params['sources'][SourceID]['name'],params['inputs']['method']))
    N = len(list(params['maps'].keys()))
    S = np.zeros(N)
    Se = np.zeros(N)
    f = np.zeros(N)
    n=0

    for MapID in tqdm(list(params['maps'].keys())):

        map_param = params['maps'][MapID]
        map_healpix = AstroMap(map_param['name'],
                               map_param['freqGHz'],1,4.5,
                               map_param['file'],
                               wcs='HPX')
        map_healpix.KelvinToJansky()
        map_wcs = map_healpix.reproject(new_wcs,new_shape)
        map_wcs.Beam2Pix()
        result = apPhotometry(map_wcs,
                              params['sources'][SourceID]['coords'],
                              method=params['inputs']['method'],
                              plotout=pdf)
        S[n],Se[n] = result.out
        f[n] = map_param['freqGHz']
        n+=1
        del map_healpix
        del result
    plt.figure()
    ax = plt.subplot()
    ax.errorbar(f,S,yerr=Se,xerr=np.ones(Se.shape)*0.5,
                 marker='x',ls='')
    ax.set_xlabel(r'$\nu$ [GHz]')
    ax.set_ylabel(r'$\log(S_\nu)$ []')
    plt.yscale('log')
    plt.suptitle('{} COMAP Source Spectrum'.format(params['sources'][SourceID]['name']))
    plt.subplots_adjust(left=0.18,right=0.95)
    pdf.savefig()

    outfile = h5py.File('Source:_{}.hd5'.format(params['sources'][SourceID]['name']),'a')
    if 'comap' not in outfile:
        outfile.create_group('comap')
    comap = outfile['comap']
    if params['inputs']['method'] in comap:
        del comap[params['inputs']['method']]
    comap.create_group(params['inputs']['method'])
    dir = comap[params['inputs']['method']]
    dir['freqGhz'] = f
    dir['Sv'] = S
    dir['Sv_err'] = Se
    outfile.close()
    pdf.close()
