;---------------------------------------------------------------
; Fill NaN nodes of a detector angle grid by a plane fit
; (z = c0 + c1*x + c2*y) using the valid nodes.
; Valid nodes are kept as they are.
;---------------------------------------------------------------
FUNCTION fill_grid_plane, grid

out = grid
dims = SIZE(grid, /DIMENSIONS)
nx = dims[0] & ny = dims[1]

xx = REBIN(FINDGEN(nx), nx, ny)               ;column index
yy = REBIN(TRANSPOSE(FINDGEN(ny)), nx, ny)    ;row index

ok  = WHERE(FINITE(grid), n_ok)
bad = WHERE(~FINITE(grid), n_bad)
IF (n_ok LT 3) OR (n_bad EQ 0) THEN RETURN, out

coef = REGRESS(TRANSPOSE([[xx[ok]], [yy[ok]]]), grid[ok], CONST=c0)
out[bad] = c0 + coef[0] * xx[bad] + coef[1] * yy[bad]

RETURN, out
END


PRO read_s2_l1c,    data_dir, sat_type, date, yyyy, mm, dd, time, hh, mi, ss, $
                    lat, lon, sza, saa, vza, vaa, raa, toa_refls, $
                    minwl, maxwl, cwl, BANDS=bands

;Bands to process (default: all 13 bands)
IF N_ELEMENTS(bands) EQ 0 THEN bands = INDGEN(13)
n_out = N_ELEMENTS(bands)

;Sentinel-2 band names (index = bandId)
band_names = ['B01','B02','B03','B04','B05','B06','B07', $
              'B08','B8A','B09','B10','B11','B12']

mtd_file = data_dir + '/MTD_MSIL1C.xml'
mtd = XML_PARSE(mtd_file)
mtd_tl_file = (FILE_SEARCH(data_dir + '/GRANULE/L1C_*/MTD_TL.xml'))[0]
mtd_tl = XML_PARSE(mtd_tl_file)

;SPACECRAFT_NAME
sat_name = mtd['n1:Level-1C_User_Product','n1:General_Info',$
               'Product_Info','Datatake','SPACECRAFT_NAME']
sat_type = STRMID(sat_name, STRLEN(sat_name)-1, 1)

;UTM Zone
utm_zone_node = mtd_tl['n1:Level-1C_Tile_ID','n1:Geometric_Info',$
                       'Tile_Geocoding','HORIZONTAL_CS_NAME']
utm_zone = (STRSPLIT(utm_zone_node, ' ', /EXTRACT))[4]

;Sensing time (UTC)
sen_time = mtd_tl['n1:Level-1C_Tile_ID','n1:General_Info',$
                  'SENSING_TIME', '#text']
date = (STRSPLIT(sen_time, 'T', /EXTRACT))[0]
yyyy = (STRSPLIT(date, '-', /EXTRACT))[0]
mm   = (STRSPLIT(date, '-', /EXTRACT))[1]
dd   = (STRSPLIT(date, '-', /EXTRACT))[2]
time = (STRSPLIT(((STRSPLIT(sen_time, 'T', /EXTRACT))[1]),'.', /EXTRACT))[0]
hh   = (STRSPLIT(time, ':', /EXTRACT))[0]
mi   = (STRSPLIT(time, ':', /EXTRACT))[1]
ss   = (STRSPLIT(time, ':', /EXTRACT))[2]

;Read lon/lat (10m res)
size_node = mtd_tl['n1:Level-1C_Tile_ID','n1:Geometric_Info',$
                   'Tile_Geocoding','Size']
nx_10 = LONG(size_node[0, 'NCOLS'])
ny_10 = LONG(size_node[0, 'NROWS'])

pos_node = mtd_tl['n1:Level-1C_Tile_ID','n1:Geometric_Info',$
                  'Tile_Geocoding', 'Geoposition']
ulx_10  = FLOAT(pos_node[0, 'ULX'])
uly_10  = FLOAT(pos_node[0, 'ULY'])
xdim_10 = FLOAT(pos_node[0, 'XDIM'])
ydim_10 = FLOAT(pos_node[0, 'YDIM'])

;Pixel center coordinates (ULX/ULY = upper-left corner of the tile)
x_utm_10 = ulx_10 + (DINDGEN(nx_10) + 0.5) * xdim_10
y_utm_10 = uly_10 + (DINDGEN(ny_10) + 0.5) * ydim_10

lon_10 = FLTARR(nx_10, ny_10) & lon_10[*] = !VALUES.F_NAN
lat_10 = FLTARR(nx_10, ny_10) & lat_10[*] = !VALUES.F_NAN

UTM2LATLON, REBIN(TEMPORARY(x_utm_10), nx_10, ny_10, /SAMPLE),$
            REBIN(TRANSPOSE(TEMPORARY(y_utm_10)), nx_10, ny_10, /SAMPLE),$
            lat_10, lon_10, ZONE=FIX(STRMID(utm_zone, 0, 2))

lat = TEMPORARY(lat_10) & lon = TEMPORARY(lon_10)

;Angle grid coordinates of 10 m pixel centers
;(angle grid: 5000 m step, starting from the tile upper-left corner)
step = 5000.0
gx = (DINDGEN(nx_10) + 0.5) * xdim_10 / step
gy = (DINDGEN(ny_10) + 0.5) * ABS(ydim_10) / step

;Read SZA
sza_node = mtd_tl['n1:Level-1C_Tile_ID','n1:Geometric_Info',$
                  'Tile_Angles',$
                  'Sun_Angles_Grid',$
                  'Zenith',$
                  'Values_List']
sza_ny = n_elements(sza_node['VALUES'])
sza_nx = n_elements(STRSPLIT((sza_node['VALUES'])[0], ' ', /EXTRACT))

sza = FLTARR(sza_nx, sza_ny) & sza[*] = !VALUES.F_NAN

FOR i = 0, sza_ny - 1 DO BEGIN
    sza[*,i] = FLOAT(STRSPLIT((sza_node['VALUES'])[i], /EXTRACT))
ENDFOR
PRINT, '# Completed to read SZA '

;Resampling SZA
sza = INTERPOLATE(sza, gx, gy, /GRID)
PRINT, '# Completed to resample SZA '

;Read SAA
saa_node = mtd_tl['n1:Level-1C_Tile_ID','n1:Geometric_Info',$
                  'Tile_Angles',$
                  'Sun_Angles_Grid',$
                  'Azimuth',$
                  'Values_List']
saa_ny = n_elements(saa_node['VALUES'])
saa_nx = n_elements(STRSPLIT((saa_node['VALUES'])[0], ' ', /EXTRACT))

saa = FLTARR(saa_nx, saa_ny) & saa[*] = !VALUES.F_NAN

FOR i = 0, saa_ny - 1 DO BEGIN
    saa[*,i] = FLOAT(STRSPLIT((saa_node['VALUES'])[i], /EXTRACT))
ENDFOR
PRINT, '# Completed to read SAA '

;Resampling SAA
saa = INTERPOLATE(saa, gx, gy, /GRID)
PRINT, '# Completed to resample SAA '

;Read VZA
vza_node = mtd_tl['n1:Level-1C_Tile_ID','n1:Geometric_Info',$
                  'Tile_Angles',$
                  'Viewing_Incidence_Angles_Grids']
vza_nx = n_elements((vza_node[0])['Zenith','Values_List','VALUES'])
vza_ny = n_elements(STRSPLIT(((vza_node[0])['Zenith','Values_List','VALUES'])[0], ' ', /EXTRACT))

n_vza_node = N_ELEMENTS(vza_node)

    ;find max,min Id
    min_b_id = 99 & max_b_id = 0
    min_d_id = 99 & max_d_id = 0

    FOR i = 0, n_vza_node - 1 DO BEGIN
        b_id = FIX((vza_node[i])['%bandId'])
        d_id = FIX((vza_node[i])['%detectorId'])

        IF b_id LT min_b_id THEN min_b_id = b_id
        IF b_id GT max_b_id THEN max_b_id = b_id
        IF d_id LT min_d_id THEN min_d_id = d_id
        IF d_id GT max_d_id THEN max_d_id = d_id
    ENDFOR

n_b_id = max_b_id - min_b_id + 1
n_d_id = max_d_id - min_d_id + 1

vza_org = FLTARR(vza_nx, vza_ny, n_b_id, n_d_id) ;nx,ny,bands,detectors
vza_org[*] = !VALUES.F_NAN

FOR i = 0, n_vza_node - 1 DO BEGIN
    b_id = FIX((vza_node[i])['%bandId'])
    d_id = FIX((vza_node[i])['%detectorId'])
    
    values_node = (vza_node[i])['Zenith','Values_List']
    FOR j = 0, vza_ny - 1 DO BEGIN
        vza_org[*,j,b_id,d_id-min_d_id] = FLOAT(STRSPLIT((values_node['VALUES'])[j], /EXTRACT))
    ENDFOR
ENDFOR
PRINT, '# Completed to read VZA '

;Read VAA
vaa_node = mtd_tl['n1:Level-1C_Tile_ID','n1:Geometric_Info',$
                  'Tile_Angles',$
                  'Viewing_Incidence_Angles_Grids']
vaa_nx = n_elements((vaa_node[0])['Azimuth','Values_List','VALUES'])
vaa_ny = n_elements(STRSPLIT(((vaa_node[0])['Azimuth','Values_List','VALUES'])[0], ' ', /EXTRACT))

n_vaa_node = N_ELEMENTS(vaa_node)

    ;find max,min Id
    min_b_id = 99 & max_b_id = 0
    min_d_id = 99 & max_d_id = 0

    FOR i = 0, n_vaa_node - 1 DO BEGIN
        b_id = FIX((vaa_node[i])['%bandId'])
        d_id = FIX((vaa_node[i])['%detectorId'])
        
        IF b_id LT min_b_id THEN min_b_id = b_id
        IF b_id GT max_b_id THEN max_b_id = b_id
        IF d_id LT min_d_id THEN min_d_id = d_id
        IF d_id GT max_d_id THEN max_d_id = d_id
    ENDFOR

n_b_id = max_b_id - min_b_id + 1
n_d_id = max_d_id - min_d_id + 1

vaa_org = FLTARR(vaa_nx, vaa_ny, n_b_id, n_d_id) ;nx,ny,bands,detectors
vaa_org[*] = !VALUES.F_NAN

FOR i = 0, n_vaa_node - 1 DO BEGIN
    b_id = FIX((vaa_node[i])['%bandId'])
    d_id = FIX((vaa_node[i])['%detectorId'])
    
    values_node = (vaa_node[i])['Azimuth','Values_List']
    FOR j = 0, vaa_ny - 1 DO BEGIN
        vaa_org[*,j,b_id,d_id-min_d_id] = FLOAT(STRSPLIT((values_node['VALUES'])[j], /EXTRACT))
    ENDFOR
ENDFOR
PRINT, '# Completed to read VAA '

;Per-pixel VZA/VAA using detector footprint masks (MSK_DETFOO)
;  - each pixel takes the angle of the detector that observed it
;  - each detector grid: NaN nodes filled by plane fit, then
;    interpolated to 10 m pixel centers
vza = FLTARR(nx_10, ny_10, n_out) & vza[*] = !VALUES.F_NAN
vaa = FLTARR(nx_10, ny_10, n_out) & vaa[*] = !VALUES.F_NAN

FOR k = 0, n_out - 1 DO BEGIN
    b = bands[k]

    ;Read detector footprint mask (pixel value = detector number, 0 = no data)
    msk_file = FILE_SEARCH(data_dir + '/GRANULE/*/QI_DATA/MSK_DETFOO_' + band_names[b] + '.jp2', COUNT=nf)
    IF nf EQ 0 THEN MESSAGE, 'Detector footprint mask not found: ' + band_names[b]

    oJP2 = OBJ_NEW('IDLffJPEG2000', msk_file[0])
    det = oJP2->GetData()
    OBJ_DESTROY, oJP2
    det = REVERSE(det, 2)

    ;20 m / 60 m mask -> 10 m grid (nearest, no interpolation)
    IF N_ELEMENTS(det[*,0]) NE nx_10 THEN det = REBIN(det, nx_10, ny_10, /SAMPLE)

    vza_band = vza[*,*,k]
    vaa_band = vaa[*,*,k]

    FOR d = 0, n_d_id - 1 DO BEGIN
        det_no = d + min_d_id
        idx = WHERE(det EQ det_no, count)
        IF count EQ 0 THEN CONTINUE

        vza_grid = fill_grid_plane(vza_org[*,*,b,d])
        vaa_grid = fill_grid_plane(vaa_org[*,*,b,d])

        vza_full = INTERPOLATE(vza_grid, gx, gy, /GRID)
        vaa_full = INTERPOLATE(vaa_grid, gx, gy, /GRID)

        vza_band[idx] = vza_full[idx]
        vaa_band[idx] = vaa_full[idx]
    ENDFOR

    vza[*,*,k] = vza_band
    vaa[*,*,k] = vaa_band
    PRINT, '# Completed to compute per-pixel VZA/VAA ' + band_names[b]
ENDFOR

!null = TEMPORARY(vza_org)
!null = TEMPORARY(vaa_org)
!null = TEMPORARY(det)

;Computing RAA
raa = FLTARR(nx_10, ny_10, n_out) ;nx,ny,bands
raa[*] = !VALUES.F_NAN

FOR i = 0, n_out - 1 DO BEGIN
    raa[*,*,i] = ABS(saa[*,*] - vaa[*,*,i])
ENDFOR
raa_idx = WHERE(raa GT 180.0, n_raa_idx)
IF n_raa_idx GT 0 THEN raa[raa_idx] = 360.0 - raa[raa_idx]
PRINT, '# Completed to compute RAA '

;QUANTIFICATION_VALUE
qa_val = mtd['n1:Level-1C_User_Product',$
             'n1:General_Info',$
             'Product_Image_Characteristics',$
             'QUANTIFICATION_VALUE',$
             '#text'] 
qa_val = FLOAT(qa_val)

;RADIO_ADD_OFFSET
offset_val_node = mtd['n1:Level-1C_User_Product',$
                   'n1:General_Info',$
                   'Product_Image_Characteristics',$
                   'Radiometric_Offset_List',$
                   'RADIO_ADD_OFFSET']
offset_vals = FLTARR(N_ELEMENTS(offset_val_node))
offset_vals[*] = !VALUES.F_NAN
FOREACH element, offset_val_node, index DO BEGIN
    offset_vals[index] = FLOAT(element['#text'])
ENDFOREACH

;SUN-EARTH Distance
d_node = mtd['n1:Level-1C_User_Product',$
             'n1:General_Info',$
             'Product_Image_Characteristics',$
             'Reflectance_Conversion',$
             'U']
d = FLOAT(d_node)

;Solar_Irradiance
F0_node = mtd['n1:Level-1C_User_Product',$
             'n1:General_Info',$
             'Product_Image_Characteristics',$
             'Reflectance_Conversion',$
             'Solar_Irradiance_List',$
             'SOLAR_IRRADIANCE']
F0s = FLTARR(N_ELEMENTS(F0_node))
F0s[*] = !VALUES.F_NAN
FOREACH element, F0_node, index DO BEGIN
    F0s[index] = FLOAT(element['#text'])
ENDFOREACH

;Spectral_Information + TOA_Reflectance
spec_info_node = mtd['n1:Level-1C_User_Product',$
                     'n1:General_Info',$
                     'Product_Image_Characteristics',$
                     'Spectral_Information_List']

minwl = FLTARR(n_out) & minwl[*] = !VALUES.F_NAN
maxwl = FLTARR(n_out) & maxwl[*] = !VALUES.F_NAN
cwl   = FLTARR(n_out) & cwl[*]   = !VALUES.F_NAN

toa_refls    = FLTARR(nx_10, ny_10, n_out)
toa_refls[*] = !VALUES.F_NAN

FOR k = 0, n_out - 1 DO BEGIN
    i = bands[k]

    ;Wavelength
    wl_node = spec_info_node['Spectral_Information', i, 'Wavelength']
    minwl[k] = FLOAT(wl_node['MIN', '#text'])
    maxwl[k] = FLOAT(wl_node['MAX', '#text'])
    cwl[k]   = FLOAT(wl_node['CENTRAL', '#text'])

    ;Find band image file
    img_data_file = FILE_SEARCH(data_dir + '/GRANULE/*/IMG_DATA/*_' + band_names[i] + '.jp2', COUNT=nf)
    IF nf EQ 0 THEN MESSAGE, 'Band file not found: ' + band_names[i]

    ;TOA Reflectance
    oJP2 = OBJ_NEW('IDLffJPEG2000', img_data_file[0])
    dn = FLOAT(oJP2->GetData())
    OBJ_DESTROY, oJP2

    ;NoData (DN = 0) -> NaN
    nodata = WHERE(dn EQ 0, n_nodata)
    IF n_nodata GT 0 THEN dn[nodata] = !VALUES.F_NAN
    PRINT, '# Completed to read ' + FILE_BASENAME(img_data_file[0])

    toa_refl = (REVERSE(dn, 2) + offset_vals[i]) / qa_val
    toa_refl = CONGRID(toa_refl, nx_10, ny_10, /INTERP, /CENTER)
    toa_refls[*,*,k] = toa_refl
    PRINT, '# Completed to compute/resample TOA Refl ' + band_names[i]

ENDFOR

END
