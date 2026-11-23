PRO map_s2_l1c_rgb
!QUIET  = 1
!EXCEPT = 0

; ==== Control Flags ====
do_Ray_cor = 'F'

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
                minwl, maxwl, cwl, BANDS=[1, 2, 3]
;bandId: B01(0), B02(1), B03(2), B04(3), B05(4), B06(5), B07(6), B08(7), B8A(8), B09(9), B10(10), B11(11), B12(12)
;Output order follows BANDS: 0 = B02, 1 = B03, 2 = B04

refl_r = REFORM(toa_refls[*, *, 2]) ;Red   (B04)
refl_g = REFORM(toa_refls[*, *, 1]) ;Green (B03)
refl_b = REFORM(toa_refls[*, *, 0]) ;Blue  (B02)

; ---- Out of range for NAN Value ----
temp = WHERE((lat LT lat_range[0]) OR (lat GT lat_range[1]) OR $
             (lon LT lon_range[0]) OR (lon GT lon_range[1]), n_temp)

IF n_temp GT 0 THEN BEGIN
    refl_r[temp] = !VALUES.F_NAN
    refl_g[temp] = !VALUES.F_NAN
    refl_b[temp] = !VALUES.F_NAN
ENDIF
    
; ==== Rayleight Correction (Bodhaine et al., 1999) ==== 
IF do_Ray_cor EQ 'T' THEN BEGIN
    PRINT, '# Start to calculate Rayleigh Optical Depth' 
        
    ;Index in toa_refls / cwl / vza / raa (BANDS=[1, 2, 3] -> 0:B02, 1:B03, 2:B04)
    ;Single scattering Rayleigh reflectance:
    ;  ref_ray = ROD * P_ray / (4 * mu0 * mu)
    ;  cos(scattering angle) = -mu0*mu - sin(sza)*sin(vza)*cos(raa)
    ;  (raa = |SAA - VAA|, raa = 0 : sun and sensor on the same side)
    ref_ray = toa_refls * 0.0
    mu0 = COS(!DTOR * sza)

    FOR k = 0, 2 DO BEGIN
        wave = cwl[k] * 0.001     ;nm -> um

        ROD = 0.0021520*(1.0455996 - 341.29061 * wave ^ (-2.0) - 0.90230850 * wave ^ 2.0) $
              / (1.0 + 0.0027059889 * wave ^ (-2.0) - 85.968563 * wave ^ 2.0)

        mu = COS(!DTOR * vza[*,*,k])
        cos_sca = -mu0 * mu - SIN(!DTOR * sza) * SIN(!DTOR * vza[*,*,k]) * COS(!DTOR * raa[*,*,k])
        P_ray = 3.0 / 4.0 * (1.0 + cos_sca ^ 2.0)

        ref_ray[*,*,k] = ROD * P_ray / (4.0 * mu0 * mu)
    ENDFOR

    refl_r = refl_r - ref_ray[*,*,2]   ;B04
    refl_g = refl_g - ref_ray[*,*,1]   ;B03
    refl_b = refl_b - ref_ray[*,*,0]   ;B02
ENDIF

; ==== mapping ====
refl_r = REVERSE(refl_r, 2)
refl_g = REVERSE(refl_g, 2)
refl_b = REVERSE(refl_b, 2)

map_name = 'S2' + sat_type + '/MSI L1C RGB' +$
           ' (' + yyyy + mm + dd + ' ' + hh + ':' + mi + ' UTC)'
   
IF do_Ray_cor EQ 'T' THEN BEGIN
    file_name = yyyy + mm + dd + '_' + hh + mi + ss +$
                '_S2' + sat_type + '_MSI_L1C_RGB_rc'
ENDIF ELSE BEGIN
    file_name = yyyy + mm + dd + '_' + hh + mi + ss +$
                '_S2' + sat_type + '_MSI_L1C_RGB'
ENDELSE
 
save_path = output_dir + 'RGB/'

PRINT, '# Completed to establish IDL-Python Bridge'
    
void    = Python.Run("import sys; sys.path.append('" + home_dir + "source')")
map_mod = Python.Import('map_rgb')
void    = map_mod.draw_rgb_map(lon, lat, Python.List(lon_range), Python.List(lat_range), $
                               refl_r, refl_g, refl_b, map_name, file_name, save_path)

PRINT, '# Completed to save : ' + file_name
        
stop
END
