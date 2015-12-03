function traitement_image,tab,tab_dark,tab_flat;,tab_flat_dark
	;Usual treatment on a raw image
	;
	result=(tab-tab_dark)/(tab_flat)
	return, result
	end

function lecture_image,filename,result
	tab=readfits(filename,header)
	tab=uint(tab)
	tab=double(tab)
	return, tab
	end

FUNCTION centre_gravite,fenetre,ecart
	; Détermine le centre de gravité d'une raie
	;Entrées : fenêtre découpée dans le FITS, Y1 ordonnée du profil, ecart au centre
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

FUNCTION pix_to_lambda,pix_raie1,pix_raie2
	;Etalonnage en longueur d'onde du spectre en utilisant les centres de gravités des raies telluriques la plus à gauche et celle juste après la deuxième raie du Fe
	;lambda_final retourne un tableau dont les composantes sont les longueurs d'onde correspondant aux pixels initiaux
	delta_pix=pix_raie2-pix_raie1
	lambda_raie1=6298.45
	lambda_raie2=6302.76
	delta_lambda=lambda_raie2-lambda_raie1
	pix_lambda=delta_lambda/delta_pix
	tot_pix=688.
	tot_lambda=tot_pix*pix_lambda
	lambda0=lambda_raie1-pix_raie1*pix_lambda
	pix_total=findgen(688)	
	lambda_final=pix_lambda*pix_total+lambda0
	RETURN,lambda_final
END

pro traitement_obs,seq_data,seq_dark,seq_flat;,seq_flat_dark
	;
	;sequence='#number of the sequences' (data,dark,flat,dark_flat)
	;the .sav file contains the array with all the processed images, the sequence numbers, the dimensions of the final array and the total number of fits files
	;example:
	;traitement_obs,'002','003','004','005'
	;restore,'p002.sav'
	;

	chemin='0928/p'+string(seq_data)+'*.fts'
	chemin_flat='0928/p'+string(seq_flat)+'*.fts'	
	chemin_dark='0928/p'+string(seq_dark)+'*.fts'

	fich1=file_search(chemin,count=nf1)
	fich2=file_search(chemin_flat,count=nf2)
	fich3=file_search(chemin_dark,count=nf3)

	dim=size(fich1)
	dimx=dim(1)	

	;
	; Traitement des images
	;
	tab_iplusv=dblarr(688,130,nf1/2.)
	tab_imoinsv=dblarr(688,130,nf1/2.)
	tab_flat=dblarr(688,520,nf2)
	tab_dark=dblarr(688,520,nf3)

	for i=0,nf2-1,1 do begin
		flat=lecture_image(fich2(i))
		tab_flat[*,*,i]=flat
	endfor
	tab_flat=total(tab_flat,3)
	tab_flat_1=tab_flat[*,125:255]
	flat_moy=total(tab_flat_1,2)
	for j=0,130 do begin
		for i=0,687 do begin
			tab_flat_1[i,j]=tab_flat_1[i,j]/flat_moy[i]
		endfor
	endfor
	help,tab_flat_1

	for i=0,nf3-1,1 do begin
		dark=lecture_image(fich3(i))
		tab_dark[*,*,i]=dark
	endfor
	tab_dark=total(tab_dark,3)
	tab_dark_1=tab_dark[*,125:255]

	for i=0,nf1-1,2 do begin
		tab=lecture_image(fich1(i))
		tab_1=tab[*,125:255]
		tab_iplusv[*,*,i/2]=traitement_image(tab_1,tab_dark_1,tab_flat_1)
	endfor
	for i=1,nf1-1,2 do begin
		tab=lecture_image(fich1(i))
		tab_1=tab[*,125:255]
		tab_imoinsv[*,*,(i-1)/2]=traitement_image(tab_1,tab_dark_1,tab_flat_1)
	endfor
	
	;
	; Calcul de V, I et V/I	et moyenne sur toutes les images
	;
	i=(tab_iplusv+tab_imoinsv)
	v=(tab_iplusv-tab_imoinsv)
	vsuri=v/i
	i=total(i,3)
	v=total(v,3)
	vsuri=total(vsuri,3)

	;
	; Sélection du milieu de la tâche pour ensuite intégrer V, I et V/I sur l'intervalle [y-5,y+5]	
	;
	window,0,xs=688,ys=130
	tvscl,i(*,*,0)
	print,"Cliquer sur la ligne la plus au milieu de l'une des deux bandes noires représentant les tâches solaires"
	cursor,x0,y0,/device
	wdelete,0
	i_fenetre=total(i(*,y0-5:y0+5),2)
	v_fenetre=total(v(*,y0-5:y0+5),2)
	vsuri_fenetre=total(vsuri(*,y0-5:y0+5),2)
	
	;
	; Sélection de la zone à étudier
	;
	window,0,xs=640,ys=480
	cgplot,vsuri_fenetre
	print,"Sélection de la zone à étudier"
	print,"Cliquer pour définir borne inférieure"
	cursor,x_inf,y1
	x_inf=uint(x_inf)
	;print,x_inf
	print,"Cliquer pour définir borne supérieure"
	cursor,x_sup,y2
	x_sup=uint(x_sup)
	;print,x_sup
	vsuri_fenetre_1=vsuri_fenetre(x_inf:x_sup)
	v_fenetre_1=v_fenetre(x_inf:x_sup)
	wdelete,0

	;
	; Tracé des profils I et V, puis sélection de la raie du fer à étudier
	;
	window,0,xs=640,ys=480
	cgplot,v_fenetre_1,xtitle="Longueur d'onde (pixels)",ytitle='Stokes V'

	window,1,xs=640,ys=480
	i_fenetre_1=i_fenetre(x_inf:x_sup)
	cgplot,i_fenetre_1,xtitle="Longueur d'onde (pixels)",ytitle='Stokes I'
	print,"Sélection de la raie à étudier"
	print,"Cliquer pour définir borne inférieure"
	cursor,x_min,y3
	x_min=uint(x_min)
	print,"Cliquer pour définir borne supérieure"
	cursor,x_max,y4
	x_max=uint(x_max)
	
	;
	; Test pour savoir quelle raie a été sélectionnée
	;
	test_raie=0
	while test_raie EQ 0 do begin
		print,format='($,"Quelle raie a été sélectionnée ?		1 -> première raie du fer					2 -> deuxième raie du fer")'
		read,test_raie
	endwhile
	print,test_raie
	wdelete,0
	wdelete,1

	;
	; Calcul des quantités nécessaires pour trouver B//
	;
	window,0,xs=640,ys=480
	cgplot,v_fenetre_1(x_min:x_max),xtitle="Longueur d'onde (pixels)",ytitle='Stokes V'
	v_max=max(v_fenetre_1(x_min:x_max))
	v_min=min(v_fenetre_1(x_min:x_max))
	;print,v_max,v_min
	index_v_max=where(v_fenetre_1(x_min:x_max) EQ v_max)
	index_v_min=where(v_fenetre_1(x_min:x_max) EQ v_min)
	delta_lambda=0.012*abs(index_v_min-index_v_max)/2. ; Calcul de delta_lambda
	;print,index_v_max,index_v_min,delta_lambda
	;print,i_fenetre_1(index_v_max+x_min),v_fenetre_1(index_v_max+x_min),vsuri_fenetre_1(index_v_max+x_min),v_fenetre_1(index_v_max+x_min)/i_fenetre_1(index_v_max+x_min)

	window,1,xs=640,ys=480
	cgplot,i_fenetre_1(x_min:x_max),xtitle="Longueur d'onde (pixels)",ytitle='Stokes I'
	print,"Sélection de l'intensité du continu"
	cursor,x1,Ic
	I_min=min(i_fenetre_1(x_min:x_max))
	delta_I=Ic-I_min ; Calcul de delta_I 
	;print,Ic,I_min,delta_I
	wdelete,0
	wdelete,1
	
	;
	; Définition des constantes nécessaires au calcul de B//
	;
	g_star_1=double(1.67)
	g_star_2=double(2.5)
	lambda1=double(6301.5)
	lambda2=double(6302.5)
	i_max=i_fenetre_1(index_v_max+x_min)
	vsuri_value=v_max/i_max

	;
	; Sélection des bonnes constantes
	;
	g_star=0.
	lambda0=0.
	if test_raie eq 1 then begin
		g_star=g_star_1
		lambda0=lambda1
	endif else begin
		g_star=g_star_2
		lambda0=lambda2
	endelse
	;print,g_star,lambda0

	;
	; Calcul de B//
	;
	r=delta_I/Ic
	A=r*exp(-1/2.)/((1-r*exp(-1/2.))*delta_lambda)
	B=4.67e-13*g_star*lambda0*lambda0
	C=A*B
	B_para=vsuri_value/C

	print, "Le champ magnétique B// vaut"+string(B_para)+" Gauss"


end
