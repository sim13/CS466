import matplotlib.pyplot as plt
import numpy as np
# match =[]

def extractor(range_start, range_end):
	match_all = []
	entropy_all = []
	entropy_rel_all = []
	for param in xrange(range_start, range_end):
		match = []
		entropy = []
		entropy_rel = []
		for i in xrange(1,11):
			path = 'Experiment Params'+ '/Experiment Params ' + str(param) + '/Dataset' + str(i)
			file_out = open(str(path) + "/eval.txt" , 'r' )
			for m ,line in enumerate(file_out):
				if m ==0:
					entropy.append(line.strip())
				if m ==1:
					match.append(line.strip())
				if m ==2:
					entropy_rel.append(line.strip())
				# Calculating percentage match
			file_out.close()
		match_all.append(match)
		entropy_all.append(entropy)
		entropy_rel_all.append(entropy_rel)
	return {"match": match_all, "rel_entropy": entropy_all, "information_content": entropy_rel_all}



def calculateMeanSD(nested_list): 
	mean_stack = []
	sd_stack = []
	for i in nested_list:
		mean_stack.append(np.mean(map(float, i)))
		sd_stack.append(np.std(map(float, i)))
	return {"mean": mean_stack, "sd": sd_stack}

def calculateStats(set):
	#_CALCS
	match_stats 				= calculateMeanSD(set["match"])
	rel_entropy_stats 			= calculateMeanSD(set["rel_entropy"])
	information_content_stats 	= calculateMeanSD(set["information_content"])
	return {"match_stats": match_stats, "rel_entropy_stats": rel_entropy_stats, "information_content_stats": information_content_stats}

def readjust(set, append_set, isML):
	# print set
	# print append_set
	if isML:
		newSet = {}
		newSet["rel_entropy"] = [set["rel_entropy"][1], set["rel_entropy"][0], append_set["rel_entropy"][0]]
		newSet["match"] = [set["match"][1], set["match"][0], append_set["match"][0]]
		newSet["information_content"] = [set["information_content"][1], set["information_content"][0], append_set["information_content"][0]]
		return newSet
	else:
		newSet = {}
		newSet["rel_entropy"] = [set["rel_entropy"][0], append_set["rel_entropy"][0], set["rel_entropy"][1]]
		newSet["match"] = [set["match"][0], append_set["match"][0], set["match"][1]]
		newSet["information_content"] = [set["information_content"][0], append_set["information_content"][0], set["information_content"][1]]
		return newSet

set1 = extractor(1,4)
set2 = extractor(4,6)
set3 = extractor(6,8)

set2_readjusted = readjust(set2, extractor(2,3), True)
set3_readjusted = readjust(set3, extractor(2,3), False)

set1_stats = calculateStats(set1)
set2_stats = calculateStats(set2_readjusted)
set3_stats = calculateStats(set3_readjusted)

print set1_stats
# print set2_stats
# print set3_stats

file_eval = open('Evaluation.txt' , 'w')
file_eval.write( 'Mean Match for NM as it goes from 0, 1,2, ' + str(set1_stats['match_stats']['mean']) +'\n' )

file_eval.write( 'Std Error Match for NM as it goes from 0, 1,2, ' + str(set1_stats['match_stats']['sd'])+'\n' )

file_eval.write( 'Mean Match for ML as it goes from 6,7,8, ' + str(set2_stats['match_stats']['mean']) +'\n')

file_eval.write( 'Std Error for ML as it goes from 6,7,8, ' + str(set2_stats['match_stats']['sd']) +'\n')

file_eval.write( 'Mean Match for SC as it goes from 5,10,20, ' + str(set3_stats['match_stats']['mean']) +'\n')

file_eval.write( 'Std Error for SC as it goes from 5,10,20, ' + str(set3_stats['match_stats']['sd']) +'\n')

file_eval.write( 'Mean rel_entropy for NM as it goes from 0, 1,2, ' + str(set1_stats['rel_entropy_stats']['mean']) +'\n')

file_eval.write( 'Std Error for NM for relative entropy as it goes from 0, 1,2, ' + str(set1_stats['rel_entropy_stats']['sd']) +'\n')

file_eval.write( 'Mean rel_entropy for  ML as it goes from 6,7,8, ' + str(set2_stats['rel_entropy_stats']['mean']) +'\n')

file_eval.write( 'Std Error rel_entropy for  ML as it goes from 6,7,8, ' + str(set2_stats['rel_entropy_stats']['sd']) +'\n')

file_eval.write( 'Mean rel_entropy for SC as it goes from 5,10,20, ' + str(set3_stats['rel_entropy_stats']['mean']) +'\n')

file_eval.write( 'Std Error rel_entropyfor SC as it goes from 5,10,20, ' + str(set3_stats['rel_entropy_stats']['sd']) +'\n')

file_eval.write( 'Mean information_content  for NM as it goes from 0, 1,2, ' + str(set1_stats['information_content_stats']['mean']) +'\n' )

file_eval.write( 'Std Error information_content for NM as it goes from 0, 1,2, ' + str(set1_stats['information_content_stats']['sd'])+'\n' )

file_eval.write( 'Mean information_content  for ML as it goes from 6,7,8, ' + str(set2_stats['information_content_stats']['mean']) +'\n')

file_eval.write( 'Std Error -information_content  for ML as it goes from 6,7,8, ' + str(set2_stats['information_content_stats']['sd']) +'\n')

file_eval.write( 'Mean information_content for SC as it goes from 5,10,20, ' + str(set3_stats['information_content_stats']['mean']) +'\n')

file_eval.write( 'Std Error of information_content for SC as it goes from 5,10,20, ' + str(set3_stats['information_content_stats']['sd']) +'\n')
# 	file_eval.write('Std Error from Mean Match Percentage for'+ str(var)+ 'is ='+ + str(snip_std) + '\n')	
file_eval.close()


# plt.plot([6,7,8] , snip_mean )
# # plt.plot([6,7,8] , snip_std , marker='o', color ='r',  linestyle ='--', label='Std Error' )
plt.errorbar([0,1,2] ,set1_stats['information_content_stats']['mean'] , set1_stats['information_content_stats']['sd'], linestyle='None' , marker='^' )
plt.axis([-1,5,-8,30])
plt.xlabel('NM Values')
plt.ylabel('information_content')
plt.title('information_content vs NM')
# # # plt.legend()
plt.show()

# # ylabel('Mean Percentage Match')



