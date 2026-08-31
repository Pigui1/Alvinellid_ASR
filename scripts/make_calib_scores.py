from subprocess import run, PIPE


ind_list = 'mut_file_sod.txt'
file_list = 'sod_consensus_7b6f2_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000_Repair.pdb'

with open(ind_list,'r') as fin:
	ll = fin.readlines()

names = []
for l in ll:
	if l[0] == '>':
		names.append(l[1:].strip())

print('>'+file_list)
cmd = './foldx_Mac/foldx_20261231 --command=Stability --pdb='+file_list
# print(cmd)
std = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
std = std.stdout.split('\n')
for i,l in enumerate(std):
	if 'FINISHING STABILITY ANALYSIS OPTION' in l:
		score = float(std[i-2].split('\t')[-1])
print(score)

for k in range(1,len(names)+1):
	wt = 'WT_'+file_list.split('.')[0]+'_'+str(k)+'.pdb'
	mm = file_list.split('.')[0]+'_'+str(k)+'.pdb'
	
	cmd = './foldx_Mac/foldx_20261231 --command=Stability --pdb='+wt
# 	print(cmd)
	std = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
	std = std.stdout.split('\n')
	for i,l in enumerate(std):
		if 'FINISHING STABILITY ANALYSIS OPTION' in l:
			score_wt = float(std[i-2].split('\t')[-1])
	
	cmd = './foldx_Mac/foldx_20261231 --command=Stability --pdb='+mm
# 	print(cmd)
	std = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
	std = std.stdout.split('\n')
	for i,l in enumerate(std):
		if 'FINISHING STABILITY ANALYSIS OPTION' in l:
			score_mm = float(std[i-2].split('\t')[-1])
	
	print('>'+names[k-1])
	print(score_mm, score_wt)
	
