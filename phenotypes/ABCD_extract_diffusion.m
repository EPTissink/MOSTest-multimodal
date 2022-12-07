# module load matlab/R2018b
# matlab -nosplash -nodesktop

cd '/path/data/'
RSI = {'N0','NF','ND'};
 
for metricnum = 1:length(RSI) 

	metric = RSI{metricnum}; 
	
	C = {'diffusion/ABCD_diffusion_PC_' metric '.mat'};  

	data = load(strjoin(C,''));                                    

	[uniqueA i j] = unique(data.subjidvec,'last');
	indexToUniques = find(ismember(1:numel(data.subjidvec),i));
	PCs = data.U(indexToUniques,:); 
	IDs = data.visitidvec(indexToUniques,:);
	out = array2table(PCs);                                
	out.Properties.RowNames = IDs;
	
	newnames = strcat(metric, '_',out.Properties.VariableNames);
	out.Properties.VariableNames = newnames;	
	
	if metricnum==1; all_metrics = out; else all_metrics = innerjoin(all_metrics,out,'Keys',{'Row'}); end	
done
 
O = {'diffusion/' metric '_PCs.txt'};
writetable(out,strjoin(O,''),'Delimiter','\t','Quotes',false,'WriteRowNames',true)
