from astropy.coordinates import SkyCoord
from astropy import units as u

import csv, re
import urllib.parse

# for nice line breaks in header strings using Astropy
# which removes all formatting characters beforehand
# just pad lines with spaces until natural linebreak
def pad_string_lines(string):
    string = string.replace('\t', '').replace("\\","").replace('  ', ' ')
    loc = 72
    last = 0
    while len(string)>loc:
        if string[loc]!=" ":
            br = loc
            spaces = ""
            while string[br]!=" " and br>last:
                br-=1
                spaces+=" "
            string = string[:br] + spaces + string[br:]
        last=loc
        loc+=72
    return string

def get_header_value(tile, value):
    if value in tile.header:
        return tile.header[value]
    else:
        return 'unknown'

def get_sexagesimal_string(position):
    sexagesimal = "%02d%02d%02.0f" % position.ra.hms+re.sub(r"([+-])\d",r"\1","%+d%02d%02d%02.0f" % position.dec.signed_dms)
    return sexagesimal

# todo incorporate source bnames into this
def get_mosaic_filename(position,radius,survey,filter=None, group_title=''):
    survey= survey.replace("PANSTARRS", "PanSTARRS")
    coords = get_sexagesimal_string(position)
    # note:the size as string already prints the units but remove space needed
    size   = str(radius).replace(" ", "")#re.sub(r"\.?0+$","","%f" % size)
    # if not self.name.replace("'", "").replace(" ", "").replace(",", "").replace("`", " ").replace('.','').replace('-','').replace('+','').isdigit():
    #     coords = self.name + "." + coords
    if not filter:
        filter=''
    if group_title=='MOSAIC':
        group_title = ''
    return f"{survey}_{str(group_title)}_J{coords}_s{size}_{filter}_mosaicked.fits"

# filename for download when NOT mosaiced and preserve as much info as possible
# replace sexigesimal location string with actual new center
def get_non_mosaic_filename(position, radius_arcmin, survey, baseurl, index, filter=None, group_title=""):
    survey= survey.replace("PANSTARRS", "PanSTARRS")
    radius = str(radius_arcmin).replace(" ", "")
    baseurl = urllib.parse.unquote(baseurl)
    basefile = baseurl.split('/')[-1].split('.fits')[0]
    new_coords = "J"+get_sexagesimal_string(position)
    if group_title=='None':
        group_title = ''
    if not filter:
        filter=''
    # if not self.name.replace("'", "").replace(" ", "").replace(",", "").replace("`", " ").replace('.','').replace('-','').replace('+','').isdigit():
    #     # if name is not just coords add it in
    #     new_coords = self.name+"." + new_coords
    if survey=='VLASS':
        old_center = re.search('J(.+?(?=\.))', basefile).group(0)
        new_base = basefile.replace(old_center, new_coords)
    else:
        new_base = survey + new_coords + "_" + basefile
    # if name is just coords use our standard coord string otherwise use their given name
    new_base = new_base.replace(survey, survey+"_"+str(group_title)+"_")
    if index!=0:
        return f"{new_base}_s{radius}_{filter}_img-{str(index+1)}.fits"
    return f"{new_base}_s{radius}_{filter}.fits"

def extractCoordfromString(position, is_name=False):
    '''
    example accepted formats:
        > RA,DEC or RA DEC in degrees
        > '00h42m30s', '+41d12m00s' or 00h42m30s, +41d12m00s
        > '00 42 30 +41 12 00'
        > '00:42.5 +41:12'
    if is_name:
        > The name of the object to get coordinates for, e.g. 'M42'
    '''
    if is_name:
        resolved = SkyCoord.from_name(position)
        return resolved
    # preformat string for consistency
    position = position.replace("'", "").replace(",", " ").replace("`", " ").replace('"','')
    if not ('h' in position or 'm' in position or ':' in position):
        # RA, DEC in degrees
        pos_list = [float(s) for s in re.findall(r'-?\d+\.?\d*', position)]
        if len(pos_list)==2: # if more than 2 it's an hour string
            return SkyCoord(pos_list[0], pos_list[1], unit=(u.deg, u.deg))

    return SkyCoord(position, unit=(u.hourangle, u.deg))

def readCoordsFromFile(csv_dictreader, max_batch=105):
    '''
    Takes a csv.DictReader object and parses coordinates from it
    Any columns matching RA, DEC header format are taken as coords
        if value nonempty
    secondarily if 'NAME' is in any headers then that value is evaluated
        as source name if value nonempty.
    Accepted variants of RA and Dec are:
    R.A.
    Right Ascension
    RA (J2000)
    R.A. (J2000)
    Right Ascension (J2000)
    RAJ2000
    DEC
    DEC.
    Declination
    DEC (J2000)
    DEC. (J2000)
    Declination (J2000)
    DecJ2000
    '''
    # REMOVE SPACES, ALL CAPS, REMOVE BRACKETS, remove decimals, REMOVE 'J2000'
    potential_RA = ['RA', 'RIGHTASCENSION']
    potential_DEC = ['DEC', 'DECLINATION']

    positions = []
    headers = csv_dictreader.fieldnames
    name_h=ra_h=dec_h=None
    errors = []
    for h in headers:
        # REMOVE SPACES, ALL CAPS, REMOVE BRACKETS, remove decimals, REMOVE 'J2000'
        trimmed = re.sub(r" ?[.()\ \[\]]", "", h.upper().replace('J2000',''))
        print(trimmed)
        if "NAME" in trimmed and not name_h:
            name_h = h
        elif trimmed in potential_RA and not ra_h:
            ra_h = h
        elif trimmed in potential_DEC and not dec_h:
            dec_h = h
    if ra_h==None or dec_h==None and name_h==None:
        raise Exception('invalid headers for coordinates or name in .CSV!')

    succ_count = 0
    line_num = 1
    for line in csv_dictreader:
        if succ_count>=max_batch:
            errors.append("max batch size is {} locations, rest were skipped".format(max_batch))
            break
        to_get = None
        is_name = False
        curr_name = ""
        # prioritize ra/dec over name but still call "input" the name for table
        if ra_h and dec_h:
            if line[ra_h]!='' and line[dec_h]!='':
                to_get = line[ra_h] +' '+ line[dec_h]
                curr_name = to_get
        if name_h:
            if line[name_h]!='':
                curr_name = line[name_h]
                if not to_get:
                    to_get = line[name_h]
                    is_name = True
        if to_get:
            try:
                # use coords for location query but report "name" as name entered no matter what for user to see
                positions.append({"name": curr_name, "position": extractCoordfromString(to_get, is_name)})
                succ_count+=1
            except Exception as e:
                errors.append(str(e))
        line_num += 1
    return positions, errors
