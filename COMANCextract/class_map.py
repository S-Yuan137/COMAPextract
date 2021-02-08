from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy_healpix import HEALPix
import healpy
import numpy as np
from astropy.coordinates import Galactic

class FileError(Exception):
    '''
    Raised when an error in the file or type of file used prevents the program
    from running successfully.
    '''
    pass

class AstroMap():
    '''
    Class to hold an instance of an astronomical map for the purposes of image
    analysis.
    '''

    def __init__(self,name,freq,freqrange,beamwidth,map,wcs=None):
        '''
        Initialise parameters:
        Frequency; FrequencyRange; Beamwidth; Map; WCS
        '''
        self.name = name
        self.freq = freq
        if type(freqrange) is int:
            self.freqrange = np.array((self.freq-abs(freqrange/2),
                                       self.freq+abs(freqrange/2)))
        else:
            self.freqrange = np.array(freqrange)
        self.beamwidth = beamwidth*u.arcmin
        if type(map) is np.ndarray:
            self.map = map
            self.wcs = wcs
        else:
            self.read_map(map,wcs)

    def read_map(self, map, wcs, hdrno=1):
        hdu = fits.open(map)
        if wcs == 'HPX':
            # PARSE AS HEALPIX
            self.wcs = {}
            self.wcs['nside'] = hdu[1].header['NSIDE']
            self.wcs['order'] = hdu[1].header['ORDERING']
            self.map = np.array(healpy.read_map(map))
            self.map[self.map<=-1.5e30] = np.nan
        else:
            # PARSE AS STANDARD GRID
            self.wcs = WCS(hdu[0].header)
            self.map = np.array(hdu[0].data)

        print(self.wcs)
        print(self.map)
        hdu.close()

    def coordinates(self):
        if type(self.wcs) is dict:
            hp = HEALPix(nside=self.wcs['nside'],
                         order=self.wcs['order'],
                         frame=Galactic())
            coords = hp.healpix_to_skycoord(np.arange(self.map.shape[0]))
            return coords
        else:
            x,y = np.meshgrid(np.arange(self.map.shape[1]),
                              np.arange(self.map.shape[0]),
                              sparse=True)
            coords = self.wcs.pixel_to_world(x,y)
            return coords

    def reproject(self,new_wcs,new_shape):
        '''
        Implement nearest-neighbour reprojection of data from one wcs to another
        wcs 'grid' image.
        CURRENTLY: only HEALPIX -> WCS supported
            work will be done to allow WCS -> WCS in future
        '''
        new_x,new_y = np.meshgrid(np.arange(new_shape[0]),
                                  np.arange(new_shape[1]),
                                  sparse=True)
        new_coords = new_wcs.pixel_to_world(new_x,new_y)
        if type(self.wcs) is dict:
            #HEALPIX projection
            hp = HEALPix(nside=self.wcs['nside'],
                         order=self.wcs['order'],
                         frame=Galactic())
            new_coords_hpx = hp.lonlat_to_healpix(new_coords.l,
                                                  new_coords.b)
            new_map = np.zeros(new_shape)
            new_map[:] = np.nan
            new_map = self.map[new_coords_hpx]
        else:
            print('Not yet supported')
        new_map_obj = AstroMap(self.name,
                               self.freq,
                               self.freqrange,
                               self.beamwidth,
                               new_map,
                               wcs=new_wcs)
        return new_map_obj

    def KelvinToJansky(self):
        beam = np.pi * ((4.5*np.pi)/(2*60*180))**2
        wavelength = 3e8/(self.freq*1e9)
        conv = (2*1.381e-23*beam)/(wavelength**2)
        conv_Jy = conv*1e26
        print('conv. to Jy/beam ({:.2f})'.format(beam))

        self.map = self.map*conv_Jy

    def Beam2Pix(self):
        beam = np.pi * np.radians(4.5/2)**2
        pix = np.radians(0.2)**2
        conv = pix/beam

        self.map = self.map*conv
