FUNCTION pix_to_lambda,pix_raie1,pix_raie2
	;Etalonnage en longueur d'onde du spectre en utilisant les centres de gravités des raies telluriques la plus à gauche et celle juste après la deuxième raie du Fe
	;lambda_final retourne un tableau dont les composantes sont les longueurs d'onde correspondant aux pixels initiaux
	delta_pix=pix_raie2-pix_raie1
	lambda_raie1=6298.45
	lambda_raie2=6302.76
	delta_lambda=lambda_raie2-lambda_raie1
	pix_lambda=delta_lambda/delta_pix
	tot_pix=680.
	tot_lambda=tot_pix*pix_lambda
	lambda0=lambda_raie1-pix_raie1*pix_lambda
	pix_total=findgen(680)	
	lambda_final=pix_lambda*pix_total+lambda0
	RETURN,lambda_final
END
