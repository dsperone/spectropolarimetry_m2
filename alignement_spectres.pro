function delta_intensity,data,pix_A,pix_G_bas,pix_A_bis,pix_G_haut
	;Fonction pour trouver le coefficient nécessaire pour aligner les continuums des 2 voies en utilisant les centres de gravités de la même raie tellurique
	pix_bas=pix_A+pix_G_bas
	pix_haut=pix_A_bis+pix_G_haut
	intensity_haut=data(int(pix_haut))
	intensity_bas=data(int(pix_bas))
	delta_intensity=intensity_haut/intensity_bas
	return,delta_intensity
	end

