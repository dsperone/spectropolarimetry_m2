FUNCTION centre_gravite,fenetre,ecart
	; Détermine le centre de gravité d'une raie
	; Entrées : fenêtre découpée dans le FITS, Y1 ordonnée du profil, ecart au centre
	; Sortie : lambda_G (centre de gravité de la raie)
   window,0
   ; Zoom sur la raie
   plot,fenetre,background=255,color=0
   print,'sélection grossière des bords de la raie pour zoomer dessus'
   print,'...bord gauche...'
   cursor,mini,Y2
   print,'...bord droit...'
   cursor,maxi,Y3
   ; Sélection du centre de la raie
   plot,fenetre(mini:maxi),background=255,color=0
   print,'...sélection du centre de la raie...'
   cursor,centre_raie,Y2
   ; Calcul des bords de la raie
   lambda_A = mini+centre_raie-ecart
   lambda_B = mini+centre_raie+ecart
   print,'lambda_A = ',lambda_A
   print,'lambda_B = ',lambda_B
   ; Intensité du continu pour calculer la fonction de poids
   print,'...intensité du continu...'
   cursor,X2,I_c
   wdelete,0
   ; Calcul du centre de gravité
   lambda = findgen(lambda_B - lambda_A + 1)
   I = fenetre(lambda_A:lambda_B)
   poids = 1.0 - I/I_c
   lambda_G = lambda_A + total(lambda*poids)/total(poids)
   print,'centre de gravité lambda G = ',lambda_G
   RETURN,lambda_G
END

FUNCTION read_fits,filename,result
	tab=readfits(filename,header)
	tab=uint(tab)
	tab=float(tab)
	RETURN, tab
END

FUNCTION pix_to_lambda,pix_raie1,pix_raie2
	;Etalonnage en longueur d'onde du spectre en utilisant les centres de gravités des raies telluriques la plus à gauche et la plus à droite du spectre
	;lambda_final retourne un tableau dont les composantes sont les longueurs d'onde correspondant aux pixels initiaux
	delta_pix=pix_raie2-pix_raie1
	lambda_raie1=6298.45
	lambda_raie2=6305.81
	delta_lambda=lambda_raie2-lambda_raie1
	pix_lambda=delta_lambda/delta_pix
	tot_pix=688.
	tot_lambda=tot_pix*pix_lambda
	lambda0=lambda_raie1-pix_raie1*pix_lambda
	pix_total=findgen(tot_pix)	
	lambda_final=pix_lambda*pix_total+lambda0
	RETURN,lambda_final
END

PRO main,seq_data,seq_dark,seq_flat,seq_dark_flat
   
   ; Transforme et formatte en chaine de caractères les num de séquence
   seq_data = string(seq_data, format='(I03)')
   seq_dark = string(seq_dark, format='(I03)')
   seq_flat = string(seq_flat, format='(I03)')
   seq_dark_flat = string(seq_dark_flat, format='(I03)')

   ;Cherche les fichiers avec les numéros de séq donnés
   chemin_data = file_search('0928/p'+seq_data+'*b1.fts',count=n_data)
   chemin_dark = file_search('0928/p'+seq_dark+'*x1.fts',count=n_dark)
   chemin_flat = file_search('0928/p'+seq_flat+'*y1.fts',count=n_flat)
   chemin_dark_flat = file_search('0928/p'+seq_dark_flat+'*x1.fts',count=n_dark_flat)

   ; Lit les fichiers fits et fait la moyenne si besoin
   data = read_fits(chemin_data(5))
   for i=6,7 do begin
      data = data + read_fits(chemin_data(i))
   endfor
   data = data / float(3)
   dark = read_fits(chemin_dark(0))
   flat = read_fits(chemin_flat(0))
   for i=1,n_flat-1 do begin
      flat = flat + read_fits(chemin_flat(i))
   endfor
   flat = flat / float(n_flat)
   dark_flat = read_fits(chemin_dark_flat(0))

   ; Détermine la distance entre les deux fenêtres de la grille
   window,0,xs=688,ys=520
   tvscl,flat
   print,'pas de grille ?'
   print,'...cheveu 1...'
   cursor,X1,Y1,/device
   ecart = 5.
   flat1 = flat(X1,Y1-ecart:Y1+ecart)
   print,'...cheveu 2...'
   cursor,X2,Y2,/device
   flat2 = flat(X1,Y2-ecart:Y2+ecart)
   min1 = min(flat1,indice_min1)
   min2 = min(flat2,indice_min2)
   pas_grille = abs(Y1-ecart+indice_min1 - (Y2-ecart+indice_min2))
   print,'pas de la grille :',pas_grille
   wdelete,0

   ; Traite le flat de son dark
   print,'corrige le dark de son flat'
   flat = flat - dark_flat
   ; Traite les données de leur dark
   print,'corrige les données du dark'
   data = data - dark

   ; Sélection de la zone d'intérêt
   print,"sélection de la zone d'intérêt"
   window,0,xs=688,ys=520
   tvscl,data
   window,1,xs=688,ys=520
   tvscl,flat
   print,"...bas de la fenêtre de la zone d'intéret..."
   cursor,X1,Y1,/device
   wdelete,0
   wdelete,1
   min1 = min(flat(X1,Y1-ecart:Y1+ecart),indice_min1)
   min2 = min(flat(X1,Y1+pas_grille-ecart:Y1+pas_grille+ecart),indice_min2)
   min3 = min(flat(X1,Y1+2*pas_grille-ecart:Y1+2*pas_grille+ecart),indice_min3)
   i1 = Y1-ecart+indice_min1
   i2 = Y1+pas_grille-ecart+indice_min2 - 1
   i3 = Y1+pas_grille-ecart+indice_min2
   i4 = Y1+2*pas_grille-ecart+indice_min3
   taille_bas = i2 - i1
   print,'contrôle :',taille_bas
   taille_haut = i4 - i3
   print,'contrôle :',taille_haut
   flat_bas = flat(*,i1:i2)
   flat_haut = flat(*,i3:i4)
   data_bas = data(*,i1:i2)
   data_haut = data(*,i3:i4)
   print,"zone d'intérêt sélectionnée"

   ; Traite le flat
   ;print,'traitement du flat : flat = flat / flat moyen'
   flat_moyen_bas = total(flat_bas,2) / taille_bas
   for j=0,taille_bas do begin
      for i=0,687 do begin
         flat_bas(i,j) = flat_bas(i,j) / flat_moyen_bas(i)
      endfor
   endfor
   flat_moyen_haut = total(flat_haut,2) / taille_haut
   for j=0,taille_haut do begin
      for i=0,687 do begin
         flat_haut(i,j) = flat_haut(i,j) / flat_moyen_haut(i)
      endfor
   endfor
   ; Affiche le flat (test)
   ;window,3
   ;tvscl,flat_haut

   ; Corrige les données du flat
   ;print,'corrige les données du flat'
   data_bas = data_bas / flat_bas
   data_haut = data_haut / flat_haut

   ; Alignement des continuums
   print,'Alignement des continuums'
   test_t = 0
   t = 0.0
   WHILE (test_t EQ 0) do begin
      ; Détermine le facteur de transmission
      print,'facteur de transmission ?'
      window,0
      tvscl,data_bas
      cursor,X1,Y1,/device
      cursor,X2,Y2,/device
      int1 = total(data_bas(X1:X2,Y1))
      int2 = total(data_haut(X1:X2,Y1))
      t = int1 / int2
      print,'t =',t
      wdelete,0
      print,'...vérification de t...'
      window,1
      plot,data_bas(*,Y1),linestyle=0,background=255,color=0
      oplot,t*data_haut(*,Y1),linestyle=0,color=2
      print,'continuums alignés ? Oui : 1, Non : 0'
      read,test_t
   endwhile
   wdelete,1
   data_haut = t * data_haut

   ; Calcul du centre de gravité d'une raie atm. pour les deux profils (décalage)
   ;print,'Calcul du décalage spectral entre les profils'
   ecart_atm = 4	; demi-largeur approximative d'une raie atmosphérique en pix
   ; Profil bas
   ;lambda_G_bas  = centre_gravite(data_bas,Y1,ecart_atm)
   ; Profil haut
   ;lambda_G_haut = centre_gravite(data_haut,Y1,ecart_atm)
   ;delta_lambda_G = abs(lambda_G_haut - lambda_G_bas)
   ;print,'décalage delta_lambda_G =',delta_lambda_G

   ; Calcul du décalage
   print,'Sélection de la raie la plus à gauche du spectre'
   pix1=centre_gravite(data_haut(*,Y1),ecart_atm)
   pix2=centre_gravite(data_bas(*,Y1),ecart_atm)
   delta_pix=abs(pix1-pix2)
   print,delta_pix

   ; Etalonnage en longueur d'onde du spectre
   print,"Etalonnage en longueur d'onde du spectre"
   print,'Sélection raie la plus à gauche'
   pix_raie1=centre_gravite(data_haut(*,Y1),ecart_atm)
   print,'Sélection raie juste à droite de la deuxième raie du Fe'
   pix_raie2=centre_gravite(data_haut(*,Y1),ecart_atm)
   lambda_final=pix_to_lambda(pix_raie1,pix_raie2)
   
   ; Interpolation linéaire pour appliquer le décalage
   n_pix=688./delta_pix
   lambda = delta_pix * findgen(n_pix)
   delta_lambda_1pix = lambda_final(1)-lambda_final(0)
   delta_lambda_delta_pix = delta_lambda_1pix*delta_pix
   lambda_final_interpol = lambda_final(0)+delta_lambda_delta_pix * findgen(n_pix)
   ;print,delta_lambda_1pix,delta_lambda_delta_pix,lambda(0:100)
   data_bas_interpol = interpolate(data_bas(*,Y1),lambda,cubic=-0.5)
   data_bas_interpol = data_bas_interpol(5:n_pix-1)
   data_haut_interpol = interpolate(data_haut(*,Y1),lambda,cubic=-0.5)
   data_haut_interpol = data_haut_interpol(0:n_pix-4)
   ;dim=size(data_bas_interpol)
   ;dim_tot=size(1)
   ;print,dim(1)
   ;cgplot,lambda_final_interpol,data_haut_interpol
   
   ; Vérification de la correction du décalage
   print,'Vérification du décalage'
   window,0
   V_sur_I = (data_haut_interpol - data_bas_interpol) / (data_haut_interpol + data_bas_interpol)
   plot,lambda_final_interpol,V_sur_I,background=255,color=0
   stop
   ; Affiche le résultat
   print,'Résultat'
   window,1
   plot,lambda_final_interpol,data_bas_interpol,background=255,color=0
   oplot,lambda_final_interpol,data_haut_interpol,color=0
   test_fe=0
   WHILE (test_t EQ 0) DO BEGIN
      print,'Passer à la sélection de la raie du Fe à étudier ? 0=Non, 1=Oui'
      read,test_fe
   ENDWHILE
   wdelete,0
   wdelete,1

   ; Sélection et étude de la raie du Fe qui nous intéresse
   print,'Sélection de la raie du Fe'
   window,0
   plot,data_bas_interpol,background=255,color=0
   oplot,data_haut_interpol,color=0
   print,'Bord gauche'
   cursor,x1,y1
   print,'Bord droit'
   cursor,x2,y2
   raie_etudiee_bas=data_bas_interpol(x1:x2)
   raie_etudiee_haut=data_haut_interpol(x1:x2)
   I=(raie_etudiee_haut-raie_etudiee_bas)/2.
   window,1
   ;cgplot,I,color=cgcolor('red')
   cgplot,raie_etudiee_bas,color=cgcolor('black'),linestyle=0
   oplot,raie_etudiee_haut,color=cgcolor('red'),linestyle=1




   stop

   ; Etalonnage en longueur d'onde en utilisant la largeur à mi-hauteur
   lim = 0
   print,'mi-hauteur ?'
   read,lim
   ; Raie atmosphère 1 = 6298.45 A
   raies_fwhm_bas = where(data_bas(5:40,Y1) LT lim)
   raies_fwhm_haut = where(t*data_haut(5:40,Y1) LT lim)
   min_raie1_bas = 5 + min(raies_fwhm_bas) + 0.5*(float(max(raies_fwhm_bas)-min(raies_fwhm_bas)))
   min_raie1_haut = 5 + min(raies_fwhm_haut) + 0.5*float((max(raies_fwhm_haut)-min(raies_fwhm_haut)))
   print,raies_fwhm_bas
   print,'raie 1'
   print,min_raie1_bas,min_raie1_haut
   ; Raie atmosphère 2 = 6299.22 A
   raies_fwhm_bas = where(data_bas(70:95,Y1) LT lim)
   raies_fwhm_haut = where(t*data_haut(70:95,Y1) LT lim)
   min_raie2_bas = 70 + min(raies_fwhm_bas) + 0.5*float((max(raies_fwhm_bas)-min(raies_fwhm_bas)))
   min_raie2_haut = 70 + min(raies_fwhm_haut) + 0.5*float((max(raies_fwhm_haut)-min(raies_fwhm_haut)))
   print,raies_fwhm_bas
   print,'raie 2'
   print,min_raie2_bas,min_raie2_haut
   ; Raie atmosphère 7 = 6305.81 A
   print,'limite raie 7 ?'
   read,lim
   raies_fwhm_bas = where(data_bas(620:650,Y1) LT lim)
   raies_fwhm_haut = where(t*data_haut(620:650,Y1) LT lim)
   min_raie7_bas = 620 + min(raies_fwhm_bas) + 0.5*float((max(raies_fwhm_bas)-min(raies_fwhm_bas)))
   min_raie7_haut = 620 + min(raies_fwhm_haut) + 0.5*float((max(raies_fwhm_haut)-min(raies_fwhm_haut)))
   print,raies_fwhm_bas
   print,'raie 7'
   print,min_raie7_bas,min_raie7_haut
   wdelete,2
   ; Correspondance longueur d'onde <-> pixel
   nb_pixel = min_raie7_bas - min_raie1_bas  
   nb_angstrom = 6305.81 - 6298.45
   print,nb_pixel
   pixel = nb_angstrom / nb_pixel
   print,' 1 pixel = ', pixel, 'A'


   ; Corrige le décalage entre les raies telluriques
   dim = size(data_bas)
   lambda_bas = fltarr(dim(1))
   lambda_haut = fltarr(dim(1))
   for i=0,dim(0)-1 do begin
      lambda_bas(i) = 6305.81 - min_raie7_bas * pixel + i * pixel
      lambda_haut(i) = 6305.81 - min_raie7_haut * pixel + i * pixel
   endfor   


   ; Trace le spectre pour chaque fenêtre de la grille
   window,1
   tvscl,data_bas
   print,'trace les profils I+V et I-V'
   window,2
   plot,data_bas(*,Y1),linestyle=0,background=255,color=0
   oplot,t*data_haut(*,Y1),linestyle=0,color=2
   window,4
   plot,((data_bas(*,Y1)-t*data_haut(*,Y1))/(data_bas(*,Y1)+t*data_haut(*,Y1))),background=255,color=0
   ;plot,(data_bas(1:680,Y1)+t*data_haut(0:679,Y1)),background=255,color=0

END
