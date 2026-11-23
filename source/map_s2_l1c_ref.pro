PRO map_s2_l1c_ref
!QUIET  = 1
!EXCEPT = 0

; ==== General Setup ====
; ---- Directories ----
home_dir   = '/home/mingyu911/S2_CalAC/'
L1C_dir    = home_dir + 'data/S2_L1C/'
output_dir = home_dir + 'output/'

; ---- Date & Location ----
yyyy = 2026 & yyyy_str = STRING(yyyy, FORMAT='(I04)')
mm   = 11   & mm_str   = STRING(mm,   FORMAT='(I02)')
dd   = 11   & dd_str   = STRING(dd,   FORMAT='(I02)')

tile_id = '52SDD'

lat_range = [34.249, 35.242]
lon_range = [127.901, 129.106]

; ==== Read Sentinel-2 L1C ====
data_dir = FILE_SEARCH(L1C_dir + 'S2A_MSIL1C_' + yyyy_str + mm_str + dd_str + '*' + tile_id + '*.SAFE', COUNT=n_data_dir)
IF (n_data_dir NE 1) THEN BEGIN
    PRINT, '!!!! Check data dir !!!!'
    STOP
ENDIF

PRINT, '# Started to read : ' + FILE_BASENAME(data_dir)

read_s2_l1c,    data_dir, sat_type, date, yyyy, mm, dd, time, hh, mi, ss, $
                lat, lon, sza, saa, vza, vaa, raa, toa_refls, $
                minwl, maxwl, cwl, BANDS=[3]
;bandId: B01(0), B02(1), B03(2), B04(3), B05(4), B06(5), B07(6), B08(7), B8A(8), B09(9), B10(10), B11(11), B12(12)
;Output order follows BANDS: 0 = B02, 1 = B03, 2 = B04

;Red (B04)
toa_refl = REFORM(toa_refls[*, *, 0])
wave     = REFORM(      cwl[*, *, 0]) 
wave_str = STRING(ROUND(wave), FORMAT='(I03)')

; ---- Out of range for NAN Value ----
temp = WHERE((lat LT lat_range[0]) OR (lat GT lat_range[1]) OR $
             (lon LT lon_range[0]) OR (lon GT lon_range[1]), n_temp)

IF n_temp GT 0 THEN BEGIN
    toa_refl[temp] = !VALUES.F_NAN
ENDIF

; ==== mapping ====
cmap       = 'Spectral_r'
map_para   = toa_refl
para_range = [0.0, 1.0]

map_name = 'S2' + sat_type + '/MSI L1C TOA Refl at ' + wave_str + ' nm' + $
           ' (' + yyyy + mm + dd + ' ' + hh + ':' + mi + ' UTC)'

file_name = yyyy + mm + dd + '_' + hh + mi + ss +$
            '_S2' + sat_type + '_MSI_L1C_TOA_Refl_' + wave_str + 'nm'
 
save_path = output_dir + 'TOA_Refl/'

PRINT, '# Completed to establish IDL-Python Bridge'
    
void    = Python.Run("import sys; sys.path.append('" + home_dir + "source')")
map_mod = Python.Import('map_utils')
void    = map_mod.draw_map(cmap, lon, lat, Python.List(lon_range), Python.List(lat_range), $
                           map_para, para_range, map_name, file_name, save_path)

PRINT, '# Completed to save : ' + file_name
        
stop
END
