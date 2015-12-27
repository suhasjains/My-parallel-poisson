awk ' BEGIN{} 
{
if(NR==1) 
	t=$3;
else if($3<t) 
	t=$3;
} 
END{print t} '  
