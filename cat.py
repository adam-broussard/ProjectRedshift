import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from sklearn.neighbors import BallTree
from astropy.coordinates import SkyCoord
from pandas import read_csv
import astropy.constants as c


def get_new_catalog():

    phot_cat = read_csv('./cosmos/cosmos.v1.3.8.cat', delim_whitespace = True, comment = '#')
    gal_cat = read_csv('./cosmos/cosmos.v1.3.6.fout', delim_whitespace = True, comment = '#')

    phot_cat_sel = phot_cat[['id', 'ra', 'dec', 'a_vector', 'b_vector', 'flux50_radius']]
    # gal_cat_sel = gal_cat[['id', 'z', 'lage', 'lmass', 'lsfr', 'lssfr']]
    gal_cat_sel = gal_cat[['id', 'z', 'lage', 'lmass']]

    new_cat = phot_cat_sel.merge(gal_cat_sel, on = 'id')

    new_cat['sfr'] = read_csv('./cosmos/cosmos.v1.3.6.sfr.v0.4.cat', delim_whitespace = True)['SFR_UVIR']

    # new_cat.lsfr = 10**new_cat.lsfr
    # new_cat.rename(columns = {'lsfr':'sfr'}, inplace = True)

    new_cat.lmass = 10**new_cat.lmass
    new_cat.rename(columns = {'lmass':'mass'}, inplace = True)

    # new_cat.lssfr = 10**new_cat.lssfr
    # new_cat.rename(columns = {'lssfr':'ssfr'}, inplace = True)

    new_cat = new_cat[(new_cat.sfr > 1) & (new_cat.sfr < 100)]

    tree = BallTree(new_cat[['ra', 'dec']])

    dists = tree.query(new_cat[['ra', 'dec']], k = 11)[0]

    new_cat['neighbor1_degrees'] = dists.T[1]
    new_cat['neighbor5_degrees'] = dists.T[5]
    new_cat['neighbor10_degrees'] = dists.T[10]




def make_ngc_cat():

    ngc_cat = read_csv('./NI2021.csv')

    sel = (ngc_cat.N == 'N') & (ngc_cat.S == 1) & np.isfinite(ngc_cat.MD) & np.isfinite(ngc_cat.z)

    # Convert RA and DEC to degrees

    ra_list = [f'{this_RH:.0f}h{this_RM:.0f}m{this_RS}s' for index, (this_RH, this_RM, this_RS) in ngc_cat[sel][['RH', 'RM', 'RS']].iterrows()]
    dec_list = [f'{this_V}{this_DG:.0f}d{this_DM:.0f}m{this_DS}s' for index, (this_V, this_DG, this_DM, this_DS) in ngc_cat[sel][['V', 'DG', 'DM', 'DS']].iterrows()]

    radec = SkyCoord(ra = ra_list, dec = dec_list)

    new_cat = ngc_cat[sel][['NI', 'A', 'SB', 'z', 'MD']]
    new_cat.A[new_cat.A.astype(str) == 'nan'] = ''

    new_cat['RA'] = radec.ra.to('deg').value
    new_cat['DEC'] = radec.dec.to('deg').value
    new_cat['VEL'] = (c.c * ((new_cat.z + 1)**2 - 1)/((new_cat.z + 1)**2 + 1)).to('km/s').value

    new_cat.rename(columns = {'NI':'NGC_ID', 'A':'EXTENSION_LETTER', 'SB':'SURFACE_BRIGHTNESS', 'MD':'COMOVING_DIST'}, inplace = True)

    new_cat.to_csv('processed_cat.csv', index = False)



def gen_group_data(ngroups = 20, galaxies_per_group = 20, min_outliers_per_group = 1):

    cat = read_csv('processed_cat.csv')
    cat['Ha_lambda'] = 6563. * (1+cat.z)
    indices = np.arange(len(cat))[cat.z > 0]

    outlier_indices = indices[((cat.VEL[indices] > 70*cat.COMOVING_DIST[indices] + 2500) | (cat.VEL[indices] < 70 * cat.COMOVING_DIST[indices] - 2000))]

    # galaxy_index_list = []

    with open('group_data.txt', 'w') as writefile:

        for this_group in range(ngroups):

            group_inds = list(np.random.choice(indices, size = galaxies_per_group - min_outliers_per_group, replace = False))
            group_inds.append(np.random.choice(outlier_indices, size = min_outliers_per_group, replace = False))
            group_inds.sort()
            # galaxy_index_list.append(group_inds.sort())

            writefile.write(repr(cat[['NGC_ID', 'COMOVING_DIST', 'Ha_lambda']].iloc[group_inds]))
            writefile.write('\n\n')




