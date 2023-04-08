#!/bin/bash


# step6.0
system="ASSEMBLY.pdb"

gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr -c ${system} -r ${system} -p topol.top -n index.ndx
gmx mdrun -v -deffnm step6.0_minimization


# Equilibration
cnt=1
cntmax=6

for ((i=${cnt}; i<${cntmax}+1; i++)); do
	if [ $i == "1" ]; then
		gmx grompp -f step6.${i}_equilibration.mdp -o step6.${i}_equilibration.tpr -c step6.$[${i}-1]_minimization.gro -r ${system} -p topol.top -n index.ndx
		gmx mdrun -v -deffnm step6.${i}_equilibration -nt 12
	else
		gmx grompp -f step6.${i}_equilibration.mdp -o step6.${i}_equilibration.tpr -c step6.$[${i}-1]_equilibration.gro -r ${system} -p topol.top -n index.ndx
		gmx mdrun -v -deffnm step6.${i}_equilibration -nt 12
	fi
done

# Production
#cnt=1
#cntmax=10
#
#for ((i=${cnt}; i<${cntmax}+1; i++)); do
#	if [ $i == "1" ]; then
#		gmx grompp -f step7_production.mdp -o step7_${i}.tpr -c step6.6_equilibration.gro -p topol.top -n index.ndx
#		gmx mdrun -v -deffnm step7_${i}
#	else
#		gmx grompp -f step7_production.mdp -o step7_${i}.tpr -c step7_$[${i}-1].gro -t step7_$[${i}-1].cpt -p topol.top -n index.ndx
#		gmx mdrun -v -deffnm step7_${i}
#	fi
#done
