awk -vCHROM="dm3.chrom.sizes" -vEXTEND=$EXTEND -vOFS='\t'
'BEGIN{while(getline>CHROM){chromSize[$1]=$2}}
 {chrom=$1;start=$2;end=$3;strand=$6;
  if(strand=="+"){
   end=start+EXTEND;
   if(end>chromSize[chrom]){
    end=chromSize[chrom]
   }
  };
  if(strand=="-"){
   start=end-EXTEND;
   if(start>1){
   start=1
   }
  };
  print chrom, start, end
 }
'
