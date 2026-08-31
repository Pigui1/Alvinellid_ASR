from subprocess import run, PIPE


folders = ['']

scores = []

for f in folders:
	comp = 0
	while 1:
		score = -1
		comp += 1
		cmd = './'+f+'foldxMacC11/foldx_20251231 --command=Stability --pdb=./'+f+'SOD_Anc1_2d5b2_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000_Repair_'+str(comp)+'.pdb'
		cmd_wt = './'+f+'foldxMacC11/foldx_20251231 --command=Stability --pdb=./'+f+'WT_SOD_Anc1_2d5b2_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000_Repair_'+str(comp)+'.pdb'
		try:
			std = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
			std = std.stdout.split('\n')
			for i,l in enumerate(std):
				if 'FINISHING STABILITY ANALYSIS OPTION' in l:
					score = float(std[i-2].split('\t')[-1])
			
			std = run(cmd_wt, stdout=PIPE, stderr=PIPE, shell=True, text=True)
			std = std.stdout.split('\n')
			for i,l in enumerate(std):
				if 'FINISHING STABILITY ANALYSIS OPTION' in l:
					score_wt = float(std[i-2].split('\t')[-1])
			scores.append([score,score_wt])
		except:
			break
			
		if score == -1:
			break

with open('scores.txt','w') as fout:
	for i in scores:
		fout.write(str(i[0])+'\t'+str(i[1])+'\n')
