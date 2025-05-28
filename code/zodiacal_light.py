import sys
import multiprocessing
import astropy.units as u
import healpy as hp
from healpy import Rotator
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
import copy
import zodipy
import pickle
import sys

DEG_TO_RAD=np.pi/180.0
clight=299792458.0

def compute_zodiacal_micron(wave_micron, model_type, NSIDE=256, date_obs="2025-01-01"):
    #get zodiacal light
    model = zodipy.Model(wave_micron * u.micron, name=model_type, extrapolate=True)
    pixels = np.arange(hp.nside2npix(NSIDE))
    lon,lat=hp.pix2ang(NSIDE, pixels, lonlat=True)
    skycoord = SkyCoord(lon, lat, unit=u.deg, frame="galactic", obstime=Time(date_obs))
    emission=model.evaluate(skycoord, nprocesses=multiprocessing.cpu_count())

    return emission, lon, lat

def compute_zodiacal(freq_GHz, model_type, NSIDE=256, date_obs="2025-01-01"):
    #get zodiacal light
    model = zodipy.Model(freq_GHz * u.GHz, name=model_type, extrapolate=True)
    pixels = np.arange(hp.nside2npix(NSIDE))
    lon,lat=hp.pix2ang(NSIDE, pixels, lonlat=True)
    skycoord = SkyCoord(lon, lat, unit=u.deg, frame="galactic", obstime=Time(date_obs))
    emission=model.evaluate(skycoord, nprocesses=multiprocessing.cpu_count())

    return emission, lon, lat

def make_zodiacal_mask(freq_GHz, model_type='planck18', eclip_deg=15.0, NSIDE=256):

    #get zodiacal light
    zd_light, lon_deg_G, lat_deg_G=compute_zodiacal(freq_GHz=freq_GHz, model_type=model_type, NSIDE=NSIDE)
    
    #convert from G to E all the coordinates/pixels
    r_G_to_E=Rotator(coord=['G','E'])
    #remember colatitute is 90-lat
    #remember if lon_lat=False, then use radians and input colatitude then longitude 
    all_theta, all_phi=r_G_to_E((90-lat_deg_G)*DEG_TO_RAD, lon_deg_G*DEG_TO_RAD)

    strip_ind=np.where((all_theta>(90-eclip_deg)*DEG_TO_RAD)&(all_theta<(90+eclip_deg)*DEG_TO_RAD))[0]
    
    #remember if lonlat=True, input longitude, latitude in degrees 
    #remember need galactic coordinates
    ec_mask=hp.pixelfunc.ang2pix(NSIDE, lon_deg_G[strip_ind], lat_deg_G[strip_ind], lonlat=True)

    #make a mask, masked region=zeros, unmasked=ones
    mask_ones=np.ones(hp.nside2npix(NSIDE))
    mask_ones[ec_mask]=0.0
    hp.write_map(f"../files/ecliptic_mask_{eclip_deg}deg.fits",mask_ones, overwrite=True)
    mask_ecliptic=hp.read_map(f"../files/ecliptic_mask_{eclip_deg}deg.fits") 
    #now some sanity checks
    
    #check if makes sense with the emission
    zd_light[ec_mask]=0.0
    plt.figure()
    hp.mollview(zd_light, unit="MJy/sr",coord=["G", "E"],title=f"Zodiacal light at {freq_GHz} GHz")
    hp.graticule()
    plt.savefig(f"../figs/{freq_GHz}_{model_type}_ec_mask_{eclip_deg}deg_E.pdf")
    
    #in galactic coord
    plt.figure()
    hp.mollview(zd_light, unit="MJy/sr",coord=["G"],title=f"Zodiacal light at {freq_GHz} GHz")
    hp.graticule()
    plt.savefig(f"../figs/{freq_GHz}_{model_type}_ec_mask_{eclip_deg}deg_G.pdf")

    #check actual mask
    plt.figure()
    hp.mollview(mask_ecliptic, coord=["G", "E"])
    hp.graticule()
    plt.savefig(f"../figs/ec_mask_{eclip_deg}deg_E.pdf")

    plt.figure()
    hp.mollview(mask_ecliptic, coord=["G"])
    hp.graticule()
    plt.savefig(f"../figs/ec_mask_{eclip_deg}deg_G.pdf")

def compute_zd(freq_GHz, model_type='planck18', date_obs="2025-01-01", mask_sky_1='', mask_sky_2='', NSIDE=256, plot_all=False):
    
    if len(mask_sky_1)==0:
        assert len(mask_sky_2)==0, "if no mask_sky_1 is provided, mask_sky_2 should also be an empty string"
    
    if 'planck' in model_type:
        zd_light_MJysr, lon_deg_G, lat_deg_G=compute_zodiacal(freq_GHz=freq_GHz, model_type=model_type, NSIDE=NSIDE, date_obs=date_obs)
    elif 'dirbe' in model_type:
        wave_micron=clight/(freq_GHz*1e9)*1e6
        zd_light_MJysr, lon_deg_G, lat_deg_G=compute_zodiacal_micron(wave_micron=wave_micron, model_type=model_type, NSIDE=NSIDE, date_obs=date_obs)
    if len(mask_sky_1)==0:
        if plot_all==True:
            #no masking
            plt.figure()
            hp.mollview(
            zd_light_MJysr,
            unit="MJy/sr",
            cmap="afmhot",coord=["G", "E"],
            title=f"Zodiacal light at {freq_GHz} GHz",
            )
            plt.savefig(f"../figs/{freq_GHz}_{model_type}_{date_obs}.pdf")
        return np.mean(zd_light_MJysr.value*1e6)

    elif len(mask_sky_1)>0:
        all_ind_1=np.arange(hp.nside2npix(NSIDE))
        #one mask
        if "HFI" in mask_sky_1:
            sky_mask_1=hp.read_map(mask_sky_1, field=3).astype(np.bool_)
            sky_1_unsmoothed=np.logical_not(sky_mask_1)
            sky_1=hp.ud_grade(sky_1_unsmoothed,NSIDE)
            sky_inds_1=all_ind_1[np.logical_not(sky_1)]
            mask_suffix_1="planck"
        elif "ecliptic" in mask_sky_1:
            sky_mask_1=hp.read_map(mask_sky_1).astype(np.bool_)
            sky_1=np.logical_not(sky_mask_1)
            sky_inds_1=all_ind_1[np.logical_not(sky_1)]
            mask_suffix_1=mask_sky_1.replace(".fits","").replace("../files/","")
    #print(np.shape(sky_inds_1))
    plot_zd_light=copy.deepcopy(zd_light_MJysr.value)
    mask_inds_1=all_ind_1[~np.logical_not(sky_1)]
    plot_zd_light[mask_inds_1]=hp.pixelfunc.UNSEEN
    
    if len(mask_sky_2)==0:
        masked_zd_light_MJysr=zd_light_MJysr[sky_inds_1]
        sky_frac=(len(sky_inds_1))/len(all_ind_1)
        if plot_all==True:
        
            #plot
            plt.figure()
            hp.mollview(plot_zd_light,
            unit="MJy/sr",
            cmap="afmhot",coord=["G", "E"],
            title=f"Zodiacal light at {freq_GHz} GHz fsky={sky_frac}",
            )
            plt.savefig(f"../figs/{freq_GHz}_{model_type}_{mask_suffix_1}_{date_obs}.pdf")

        return np.mean(masked_zd_light_MJysr.value*1e6)
    else:
        all_ind_2=np.arange(hp.nside2npix(NSIDE))
        if "HFI" in mask_sky_2:
            sky_mask_2=hp.read_map(mask_sky_2, field=3).astype(np.bool_)
            sky_2_unsmoothed=np.logical_not(sky_mask_2)
            sky_2=hp.ud_grade(sky_2_unsmoothed,NSIDE)
            sky_inds_2=all_ind_2[np.logical_not(sky_2)]
            mask_suffix_2="planck"
        elif "ecliptic" in mask_sky_2:
            sky_mask_2=hp.read_map(mask_sky_2).astype(np.bool_)
            sky_2=np.logical_not(sky_mask_2)
            sky_inds_2=all_ind_2[np.logical_not(sky_2)]
            mask_suffix_2=mask_sky_2.replace(".fits","").replace("../files/","")

        sky_inds_both=np.unique(np.intersect1d(sky_inds_1, sky_inds_2))
        #print(np.shape(sky_inds_both))
        masked_zd_light_MJysr=zd_light_MJysr[sky_inds_both] 
        
        #for plot
        mask_inds_2=all_ind_2[~np.logical_not(sky_2)]
        plot_zd_light[mask_inds_2]=hp.pixelfunc.UNSEEN
        
        mask_inds_both=np.unique(np.concatenate((mask_inds_1, mask_inds_2)))
        plot_zd_light_test=copy.deepcopy(zd_light_MJysr.value)
        plot_zd_light_test[mask_inds_both]=hp.pixelfunc.UNSEEN

        test_pixels_1=np.where(plot_zd_light_test!=hp.pixelfunc.UNSEEN)[0]
        test_pixels_2=np.where(plot_zd_light!=hp.pixelfunc.UNSEEN)[0]
        #print(len(test_pixels_1))
        #print(len(test_pixels_2))

        assert len(test_pixels_1)==len(test_pixels_2)
        assert np.array_equal(plot_zd_light_test,plot_zd_light)

        sky_frac=len(sky_inds_both)/len(all_ind_1)
        if plot_all==True:
            #plot
            plt.figure()
            hp.mollview(plot_zd_light,
            unit="MJy/sr",
            cmap="afmhot",coord=["G", "E"],
            title=f"Zodiacal light at {freq_GHz} GHz fsky={sky_frac}",
            )
            plt.savefig(f"../figs/{freq_GHz}_{model_type}_{mask_suffix_1}_{mask_suffix_2}_{date_obs}.pdf")

        return np.mean(masked_zd_light_MJysr.value*1e6)


