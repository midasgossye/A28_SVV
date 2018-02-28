##################################################
#### VERIFICATION CHECKER PROG ###################
##################################################


#==========================================MAIN PROG==========================================#

print '#### VERIFICATION CHECKER ####'
print
print
print '1. shear forces'
print '2. shear flows'
print '3. deflections LE/TE'
print '4. geometrical properties'
print

start = raw_input('Enter an option from the given menu:  ')

if start == '1':                                                                                        # COMPARISON OF SHEAR FORCE VALUES IN THE ANALYTICAL
                                                                                                        # AND THE NUMERICAL MODELS
    shear_f_az = float(input('Enter the shear force from the analytical model [N] in (+)z-direction:   '))        

    shear_f_nz = float(input('Enter the shear force from the numerical model [N] in (+)z-direction:   '))

    shear_f_ay = float(input('Enter the shear force from the analytical model [N] in (+)y-direction:   '))

    shear_f_ny = float(input('Enter the shear force from the numerical model [N] in (+)y-direction:   '))

    error_shear_f_z = abs((shear_f_az - shear_f_nz)/shear_f_az)*100

    error_shear_f_y = abs((shear_f_ay - shear_f_ny)/shear_f_ay)*100
    
    print 'Error for shear force in z-dir is', error_shear_f_z,'%'
    print 'Error for shear force in y-dir is', error_shear_f_y,'%'

elif start == '2':                                                                                      # COMPARISON OF SHEAR FLOW VALUES IN THE ANALYTICAL
                                                                                                        # AND THE NUMERICAL MODELS
    shear_sf_az = float(input('Enter the shear flow from the analytical model [N/mm] in (+)z-direction:   '))

    shear_sf_nz = float(input('Enter the shear flow from the numerical model [N/mm] in (+)z-direction:   '))

    shear_sf_ay = float(input('Enter the shear flow from the analytical model [N/mm] in (+)y-direction:   '))

    shear_sf_ny = float(input('Enter the shear flow from the numerical model [N/mm] in (+)y-direction:   '))

    error_shear_sf_z = abs((shear_az - shear_nz)/shear_az)*100

    error_shear_sf_y = abs((shear_ay - shear_ny)/shear_ay)*100
    
    print 'Error for shear flow in z-dir is', error_shear_sf_z,'%'
    print 'Error for shear flow in y-dir is', error_shear_sf_y,'%'

elif start == '3':                                                                                      # COMPARISON OF DEFLECTION VALUES IN THE ANALYTICAL
                                                                                                        # AND THE NUMERICAL MODELS @ LEADING AND TRAILING EDGE                                 
    deflec_a_le = float(input('Enter the delfection of the LE of the analytical model in (+)y-dir.[mm]:   '))     

    deflec_n_le = float(input('Enter the delfection of the LE of the numerical model in (+)y-dir.[mm]:   '))

    deflec_a_te = float(input('Enter the delfection of the TE of the analytical model in (+)y-dir.[mm]:   '))

    deflec_n_te = float(input('Enter the delfection of the TE of the numerical model in (+)y-dir.[mm]:   '))

    error_deflec_le = abs((deflec_a_le - deflec_n_le)/deflec_a_le)*100

    error_deflec_te = abs((deflec_a_te - deflec_n_te)/deflec_a_te)*100
    
    print 'Error for deflection in y-dir @ LE is', error_deflec_le,'%'
    print 'Error for deflection in y-dir @ TE is', error_deflec_te,'%'

elif start == '4':                                                                                      # COMPARISON OF VARIOUS GEOMETRICAL PROPERTIES OF THE 
                                                                                                        # ANALTYICAL AND NUMERICAL MODELS
    print
    print '1. centroid location'
    print '2. cross sectional area'
    print '3. moment of inertia y-dir.'
    print '4. moment of inertia z-dir.'
    print '5. moment of inertia yz-dir.'
    print
    
    start2 = raw_input('Enter any of the options given:   ')

    if start2 == '1':

        c_a = float(input('Enter the centroid location seen from (+)z-dir. @ TE for the analytical model [mm]:   '))

        c_n = float(input('Enter the centroid location seen from (+)z-dir. @ TE for the numerical model [mm]:   '))

        error_c = abs((c_a-c_n)/c_a)*100
    
        print 'Error for centroid in z-dir is', error_c,'%'

    elif start2 == '2':

        ar_a = float(input('Enter the cross sectional area for the analytical model [mm^2]:   '))

        ar_n = float(input('Enter the cross sectional area for the numerical model [mm^2]:   '))

        error_ar = abs((ar_a-ar_n)/ar_a)*100
    
        print 'Error for cross sectional area is', error_ar,'%'

    elif start2 == '3':

        moiy_a = float(input('Enter the MOI y-dir. for the analytical model [mm^4]:   '))

        moiy_n = float(input('Enter the MOI y-dir. for the numerical model [mm^4]:   '))

        error_moiy = abs((moiy_a-moiy_n)/moiy_a)*100
    
        print 'Error for MOI y-dir. is', error_moiy,'%'

    elif start2 == '4':

        moiz_a = float(input('Enter the MOI z-dir. for the analytical model [mm^4]:   '))

        moiz_n = float(input('Enter the MOI z-dir. for the numerical model [mm^4]:   '))

        error_moiz = abs((moiz_a-moiz_n)/moiz_a)*100
    
        print 'Error for MOI z-dir. is', error_moiz,'%'

    elif start2 == '5':

        moiyz_a = float(input('Enter the MOI yz-dir. for the analytical model [mm^4]:   '))

        moiyz_n = float(input('Enter the MOI yz-dir. for the numerical model [mm^4]:   '))

        error_moiyz = abs((moiyz_a-moiyz_n)/moiyz_a)*100
    
        print 'Error for MOI yz-dir. is', error_moiyz,'%'

else:
    print 'BZZT ERROR!  TRY RUNNING PROG. AGAIN'

        
    

        

    

    










    

