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

class grizy_filters(Enum):
    g = 1 
    r = 2
    i = 3
    z = 4
    y = 5


#
#   H E A D E R   F I L T E R I N G / F O R M A T I N G
#

def get_header_pretty_string(header):
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


class HeaderFilter:
    def __init__(self, header, is_add_wcs=False):
        # set the header layout order, for saved_keys, and default comments, if any.
        self.header_layout = {
            'TOP_KEYS': [
                ['SIMPLE'],
                ['BITPIX'],
                ['NAXIS',{'COMMENT': 'Number of array dimensions'}],
                [{'NAXIS_BLOCK': [
                    ['NAXIS1'],
                ]}],
                ['WCSAXES',{'VALUE': 'FK5', 'COMMENT': 'Number of WCS axes'}],
                [{'WCSAXES_BLOCK': [
                    ['CTYPE1'],
                    ['CUNIT1'],
                    ['CRVAL1'],
                    ['CRPIX1',{'COMMENT': 'Axis %d reference pixel'}],
                    ['CDELT1'],
                ]}],
                ['LATPOLE'],
                ['LONPOLE'],
                ['RADESYS'],
                ['SURVEY'],
                ['BAND']
            ],
            'BOTTOM_KEYS': [
                ['EPOCH'],
                ['DATE-OBS'],
                ['COMMENT']
            ]
        }

        self.saved_keys = list()
        self.updates = dict()
        self.header = header
        if is_add_wcs:
            self.update(WCS(header).to_header())

        # special header keys
        self.reserved_keys = ['SIMPLE', 'BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND']
        self.comment_block = ['COMMENT']



    def __set_layout(self):
        def get_keys(blocks):
            naxis   = self.header['NAXIS'] if 'NAXIS' in self.header else None
            wcsaxes = self.header['WCSAXES'] if 'WCSAXES' in self.header else None
            block_keys = list()
            for block in blocks:
                if isinstance(block[0],str):
                    block_keys.append(block[0])
                elif isinstance(block[0],dict):
                    if 'NAXIS_BLOCK' in block[0] and naxis:
                        fields = block[0]['NAXIS_BLOCK']
                        for n in range(1,naxis+1):
                            for field in fields:
                                block_keys.append(re.sub(r"\d$",f"{n}",field[0]))
                    elif 'WCSAXES_BLOCK' in block[0] and wcsaxes:
                        fields = block[0]['WCSAXES_BLOCK']
                        for n in range(1,wcsaxes+1):
                            for field in fields:
                                block_keys.append(re.sub(r"\d$",f"{n}",field[0]))
            return block_keys
        self.saved_keys = list(set(self.saved_keys))
        top_keys    = get_keys(self.header_layout['TOP_KEYS'])
        bottom_keys = get_keys(self.header_layout['BOTTOM_KEYS'])
        if len(list(set(top_keys) & set(bottom_keys))) > 0:
            # TODO: Make into warning...
            print("Whoops!")
            return self.saved_keys
        top_set    = list(set(top_keys) & set(self.saved_keys))
        bottom_set = list(set(bottom_keys) & set(self.saved_keys))
        middle_set = [k for k in self.saved_keys if not k in top_set and not k in bottom_set]

        unordered_keys = self.saved_keys.copy()
        ordered_keys = list()
        for k in top_keys:
            if k in unordered_keys:
                ordered_keys.append(k)
                unordered_keys.remove(k)
        for k in unordered_keys:
            if k in middle_set:
                ordered_keys.append(k)
        for k in middle_set:
            if k in unordered_keys:
                unordered_keys.remove(k)
        for k in bottom_keys:
            if k in unordered_keys:
                ordered_keys.append(k)
                unordered_keys.remove(k)

        if len(unordered_keys) == 0:
            self.saved_keys = ordered_keys.copy()
        else:
            print(f"{type(self).__name__}: WARNING: unable to order keys: ramainder {unordered_keys}.")

        return self.saved_keys

    def do_it(self):
        self.__set_layout()

    def update(self,updates,is_overwrite_existing=True):
        if not (updates is None):
            self.saved_keys.extend([k for k in updates.keys()])
            if is_overwrite_existing:
                self.updates.update(updates)
                self.header.update(updates)
            else:
                for k in updates.keys():
                    if k in self.header:
                        self.updates[k] = (self.header[k], self.header.comments[k])
                    else:
                        self.updates[k] = updates[k]
                        self.header[k]  = updates[k]
            self.__set_layout()
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
            self.__set_layout()
        return self.get_saved_keys()

    def get_saved_keys(self):
        return self.saved_keys

    def get_updates(self):
        return self.updates

    def get_header(self):
        return self.header
