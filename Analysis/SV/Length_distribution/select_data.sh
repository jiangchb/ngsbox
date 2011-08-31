mysql -N -u korbinian -pstephan -h ume.eb.local ath_pe -e "select length from sv_indel where type = 'PE' and sample = 'Bur-0' and support/members > 0.9" > sv_indel.pe.bur-0.txt
mysql -N -u korbinian -pstephan -h ume.eb.local ath_pe -e "select length from sv_indel where type = 'MP' and sample = 'Bur-0' and support/members > 0.9" > sv_indel.mp.bur-0.txt
mysql -N -u korbinian -pstephan -h ume.eb.local ath_pe -e "select length from sv_indel where type = 'PE' and sample = 'Col-0' and support/members > 0.9" > sv_indel.pe.col-0.txt
