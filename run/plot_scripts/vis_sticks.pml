remove solvent
select n. zn+
remove sele
select n. na+
remove sele
select DNA, resn DA+DC+DG+DT+PSD
select INT, interface
select PRO, not DNA or INT
hide cartoon, PRO
show surface, PRO
select PRO_INT, br. PRO within 5 of DNA
show sticks, PRO_INT
set cartoon_ring_mode, 3
set cartoon_nucleic_acid_color, wheat
select AT, resn DA+DT
set cartoon_ring_color, palegreen, AT
set cartoon_ring_color, lightblue, not AT
split_chain
set transparency, 0.7
set surface_color, grey90
set cartoon_color, grey90
set stick_color, grey90
set specular, 0

