import re
from enum import Enum
from astropy.wcs import WCS


#  S U V E R Y   F I L T E R S
#
# WARNING: Don't change the enum number assignements in this file
#


class wise_filters(Enum):
    w1 = 1
    w2 = 2
    w3 = 3
    w4 = 4
    # w1 = "W1 (3.4 %CE%BCm)"
    # w2 = "W2 (4.6μm)"
    # w3 = "W3 (12μm)"
    # w4 = "W4 (22μm)"


class grizy_filters(Enum):
    g = 1
    r = 2
    i = 3
    z = 4
    y = 5

class ugriz_filters(Enum):
    u = 1
    g = 2
    r = 3
    i = 4
    z = 5

class gleam_frequency(Enum):
    f1 = "072-103"
    f2 = "103-134"
    f3 = "139-170"
    f4 = "170-231"

class vlass_epoch(Enum):
    e11 = "1.1"
    e12 = "1.2"
    e21 = "2.1"

#
#   H E A D E R   F I L T E R I N G / F O R M A T I N G
#

def get_header_pretty_string(header):
    """This is routine for providing nice formatted header strings
       for pretty printing."""
    def get_value(h,k):
        return h[k][0] if isinstance(header,dict) else h[k]
    def get_comments(h,k):
        return h[k][1] if isinstance(header,dict) else h.comments[k]
    pr_str = ""
    for k in header.keys():
        hdr_value = get_value(header,k)
        if isinstance(hdr_value, str):
             hdr_vstr = "'%s'" % hdr_value
        elif isinstance(hdr_value,int):
             if k == 'SIMPLE' or k == 'EXTEND':
                 hdr_vstr = 'T' if hdr_value > 0 else 'F'
             else:
                 hdr_vstr = "%d" % hdr_value
        elif isinstance(hdr_value,float):
             hdr_vstr = "%f" % hdr_value
        else:
             hdr_vstr = f"{hdr_value}"
        pr_str += f"{k:<9} = {hdr_vstr:<21} / {get_comments(header,k)}"+"\n"
    return pr_str


def sanitize_fits_date_fields(date_obs_value):
    """This routine attempts to ensure the fits header date field conforms to
       post-Y2K stanards: i.e.,

          KEYWORD:   DATE-OBS
          REFERENCE: FITS Standard
          STATUS:    reserved
          HDU:       any
          VALUE:     string
          COMMENT:   date of the observation
          DEFINITION: The date of the observation, in the format specified in the
          FITS Standard.  The old date format was 'yy/mm/dd' and may be used only
          for dates from 1900 through 1999.  The new Y2K compliant date format is
          'yyyy-mm-dd' or 'yyyy-mm-ddTHH:MM:SS[.sss]'.

       as per, https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html.
    """
    # attempts to fix y2k issues, re.,
    # https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html,
    # and header flaws...
    date_obs = re.sub(r"\s+","",date_obs_value)
    reduced_date_obs = re.sub(r"(\s|-|/)","",date_obs_value)
    if len(reduced_date_obs) == 6: # pre-y2k
        # nb: this is not perfect: i.e., it has not been completely rationalized out.
        try:
            year_xfix     = int(reduced_date_obs[0:2])
            year_or_month = int(reduced_date_obs[2:4])
            month_or_day  = int(reduced_date_obs[4:6])
            if not (year_xfix == 19 or year_xfix == 20) and \
               (1 <= year_or_month and year_or_month <= 12) and \
               (1 <= month_or_day  and month_or_day  <= 31):
                date_obs = '19' + reduced_date_obs
            else: # ok, not y2k, it's yyyymm -- violates standards
                date_obs = reduced_date_obs + "15"
            date_obs = f"{int(date_obs[0:4])}-{int(date_obs[4:6]):02d}-{int(date_obs[6:8]):02d}"
        except:
            pass
    elif len(reduced_date_obs) == 8: # yyyymmdd not post-y2k standard
        date_obs = f"{int(reduced_date_obs[0:4])}-{int(reduced_date_obs[4:6]):02d}-{int(reduced_date_obs[6:8]):02d}"
    return date_obs


class HeaderFilter:
    """This class is used as a header filter tool for modifying, adding, formating, and extracting header infomation
       intended for the creation of new headers.

       method:
          - update(updates,is_overwrite_existing=True) - updates internal header and adds field to save list
          - save_keys(keys)  - adds header field to save list
          - get_saved_keys() - returns save list ordered according to self.header_layout
          - get_header()     - returns modified copy of header

    """
    def __init__(self, header, is_add_wcs=False):

        # set the header layout order, for saved_keys, and default comments, if any.
        # Notes:
        #    [1] Not in this list are placed between TOP_KEYS and BOTTOM_KEYS.
        #    [2] The default (i.e., if none) COMMENT fields have not been implemented.
        # TODO (Issue #6): Determine if default COMMENT are still required; if yes, implement.
        # TODO (Issue #6): Remove default VALUE field -- only second guessing... so...
        # self.header_layout = {
        #     'TOP_KEYS': [
        #         ['SIMPLE'],
        #         ['BITPIX'],
        #         ['NAXIS',{'COMMENT': 'Number of array dimensions'}],
        #         [{'NAXIS_BLOCK': [
        #             ['NAXIS1'],
        #         ]}],
        #         ['WCSAXES',{'VALUE': '???', 'COMMENT': 'Number of WCS axes'}],
        #         ['RADESYS'],
        #         ['PC1_1'],
        #         ['PC1_2'],
        #         ['PC2_1'],
        #         ['PC2_2'],
        #         [{'WCSAXES_BLOCK': [
        #             ['CTYPE1'],
        #             ['CUNIT1'],
        #             ['CRVAL1'],
        #             ['CRPIX1',{'COMMENT': 'Axis %d reference pixel'}],
        #             ['CDELT1'],
        #         ]}],
        #         ['LATPOLE'],
        #         ['LONPOLE'],
        #         ['SURVEY'],
        #         ['BAND'],
        #         ['STOKES'],
        #         ['STK_TYPE',{'COMMENT': 'Survey image stack type'}],
        #         ['STK_ID',  {'COMMENT': 'Survey image sksy cell'}],
        #         ['SKYCELL', {'COMMENT': 'Survey image sky cell'}],
        #         ['TESS_ID', {'COMMENT': 'Survey tesselation'}],
        #         ['BUNIT',{'VALUE': '???', 'COMMENT': 'Pixel flux unit'}],
        #         ['BMAJ', {'COMMENT': 'Beam major axis [deg]'}],
        #         ['BMIN', {'COMMENT': 'Beam minor axis [deg]'}],
        #         ['BPA',  {'COMMENT': 'Beam position angle'}],
        #     ],
        #     'BOTTOM_KEYS': [
        #         ['EQUINOX'],
        #         ['EPOCH'],
        #         ['MJD'],
        #         ['MJD-OBS'],
        #         ['DATE-OBS'],
        #         ['FNAME', {'COMMENT': 'Survey coadded image'}],
        #         ['IMFILE',{'COMMENT': 'ATLAS image identifier'}],
        #         ['COMMENT']
        #     ]
        # }

        self.saved_keys = list()
        self.updates = dict()
        self.header = header.copy()
        # special header keys
        self.reserved_keys = ['SIMPLE', 'BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND']
        self.comment_block = ['COMMENT']
        # sanitize the date-obs string
        if 'DATE-OBS' in header:
            self.header['DATE-OBS'] = (sanitize_fits_date_fields(self.header['DATE-OBS']), self.header.comments['DATE-OBS'])

        # Use WCS to do caluclations and extract required fields for 2D (WCSAXES) fits images, only.
        if is_add_wcs:
            wcs_header = WCS(header,naxis=2).to_header()
            keep = [
                'WCSAXES',
                'RADESYS',
                'PC1_1',
                'PC1_2',
                'PC2_1',
                'PC2_2',
                'CTYPE1',
                'CUNIT1',
                'CRVAL1',
                'CRPIX1',
                'CDELT1',
                'CTYPE2',
                'CUNIT2',
                'CRVAL2',
                'CRPIX2',
                'CDELT2',
                'LATPOLE',
                'LONPOLE',
                'RADESYS',
                'EQUINOX',
                'MJD-OBS',
                'DATE-OBS'
            ]
            for field in wcs_header:
                # if field in keep:
                #     self.update({field: (wcs_header[field], wcs_header.comments[field])})
                # TODO False means DONT OVERWRITE.... should we???
                self.update({field: (wcs_header[field], wcs_header.comments[field])}, True)
                keep.append(field)
            self.update({'WCSAXES': (2, wcs_header.comments['WCSAXES'])})
            self.save_keys(keep)




    # def __set_layout(self):
    #     """This is a private routine which is (should be) called each time self.saved_keys is updated,
    #        in order to insure it's order according the pattern defined in self.header_layout."""
    #     # helper function to extract order-keys (header fields) from the self.header_layout
    #     # top and bottom blocks.
    #     def get_keys(blocks):
    #         naxis   = self.header['NAXIS'] if 'NAXIS' in self.header else None
    #         wcsaxes = self.header['WCSAXES'] if 'WCSAXES' in self.header else None
    #         block_keys = list()
    #         for block in blocks:
    #             if isinstance(block[0],str):
    #                 block_keys.append(block[0])
    #             elif isinstance(block[0],dict):
    #                 if 'NAXIS_BLOCK' in block[0] and naxis:
    #                     fields = block[0]['NAXIS_BLOCK']
    #                     for n in range(1,naxis+1):
    #                         for field in fields:
    #                             block_keys.append(re.sub(r"\d$",f"{n}",field[0]))
    #                 elif 'WCSAXES_BLOCK' in block[0] and wcsaxes:
    #                     fields = block[0]['WCSAXES_BLOCK']
    #                     for n in range(1,wcsaxes+1):
    #                         for field in fields:
    #                             block_keys.append(re.sub(r"\d$",f"{n}",field[0]))
    #         return block_keys
        #
        # # break up the keys into top, bottom, and middle sets.
        # self.saved_keys = list(set(self.saved_keys))
        # top_keys    = get_keys(self.header_layout['TOP_KEYS'])
        # bottom_keys = get_keys(self.header_layout['BOTTOM_KEYS'])
        # if len(list(set(top_keys) & set(bottom_keys))) > 0:
        #     # TODO (Issue #6): Make into warning...
        #     print("Whoops!")
        #     return self.saved_keys
        # top_set    = list(set(top_keys) & set(self.saved_keys))
        # bottom_set = list(set(bottom_keys) & set(self.saved_keys))
        # middle_set = [k for k in self.saved_keys if not k in top_set and not k in bottom_set]
        #
        # # now lets order the self.saved_keys according to self.header_layout
        # unordered_keys = self.saved_keys.copy() # nb: must copy in place
        # ordered_keys = list()
        # for k in top_keys:
        #     if k in unordered_keys:
        #         ordered_keys.append(k)
        #         unordered_keys.remove(k)
        # for k in unordered_keys:
        #     if k in middle_set:
        #         ordered_keys.append(k)
        # for k in middle_set:
        #     if k in unordered_keys:
        #         unordered_keys.remove(k)
        # for k in bottom_keys:
        #     if k in unordered_keys:
        #         ordered_keys.append(k)
        #         unordered_keys.remove(k)
        #
        # # copy sorted keys if all accounted for
        # if len(unordered_keys) == 0:
        #     self.saved_keys = ordered_keys.copy()
        # else:
        #     print(f"{type(self).__name__}: WARNING: unable to order keys: ramainder {unordered_keys}.")
        #
        # return self.saved_keys


    def update(self,updates,is_overwrite_existing=True):
        """This is use to update the running self.header while keeping tabs of what was updated
           for later extraction, via, self.saved keys (along with self.updates that keep a record)."""
        # update if not empty
        if not (updates is None):
            self.saved_keys.extend([k for k in updates.keys()])
            if is_overwrite_existing: # upate everything
                self.updates.update(updates)
                self.header.update(updates)
            else: # don't overwrite self.header, but update if field doesn't exist.
                for k in updates.keys():
                    if k in self.header:
                        self.updates[k] = (self.header[k], self.header.comments[k])
                    else:
                        self.updates[k] = updates[k]
                        self.header[k]  = updates[k]
            # reorder keys to conform to self.header_layout
            # self.__set_layout()
        return self


    def save_keys(self,keys):
        if not (keys is None):
            if isinstance(keys,str):
                keys = [keys]
            for key in keys:
                k = key.upper()
                if not (k in self.saved_keys) and k in self.header:
                    self.updates[k] = (self.header[k], self.header.comments[k])
                    self.saved_keys.append(k)
        #     self.__set_layout()
        # return self.get_saved_keys()
        return self.saved_keys

        # FALON COMMENTED TGUS OUT RECENTLY
    # def get_saved_keys(self):
    #     """The will return a sorted list of saved header fields, as define by self.header_layout."""
    #     return self.saved_keys


    ## TODO (Issue #6): *** DEPRECATED ***
    ##       Appears to be of no value... consider removing...
    #def get_updates(self):
    #    """The will return a dict of updates."""
    #    return self.updates


    def get_header(self):
        """This returns the modified header, which can be used in conjucting with get_save_keys(),
           to provide a nicely sorted header, in accordance with self.header_layout."""
        return self.header
