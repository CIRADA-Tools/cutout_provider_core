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
        try:
           return h[k][0]
        except:
           pass
        return h[k]
    def get_comments(h,k):
        try:
            return h.comments[k]
        except:
            return h[k][1]
    pr_str = ""
    for k in header.keys():
        hdr_value = get_value(header,k)
        if isinstance(hdr_value, str):
             hdr_vstr = "'%s'" % hdr_value
        elif isinstance(hdr_value,int):
             hdr_vstr = "%d" % hdr_value
        elif isinstance(hdr_value,float):
             hdr_vstr = "%f" % hdr_value
        else:
             hdr_vstr = f"{hdr_value}"
        pr_str += f"{k:<9} = {hdr_vstr:<21} / {get_comments(header,k)}"+"\n"
    return pr_str

class HeaderFilter:
    def __init__(self,header,is_add_wcs=False):
        self.saved_keys = list()
        self.updates = dict()
        self.header = header
        if is_add_wcs:
            self.update(WCS(header).to_header())

    def update(self,updates,is_overwrite_existing=True):
        if not (updates is None):
            self.saved_keys.extend([k for k in updates.keys()])
            self.saved_keys = list(set(self.saved_keys))
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
        return self

    def get_saved_keys(self):
        return self.saved_keys

    def get_updates(self):
        return self.updates

    def get_header(self):
        return self.header
