{
    if(/Memcheck/) file=new;
    new=$0;
    if(/HEAP SUMMARY/){
	print file;
    }
    if(/definitely lost/){
	if($4!=0) print $0;
    }
    if(/still reachable/){
	if($4!=0) print $0;
    }
    if(/ERROR SUMMARY/){
	if($4!=0){
	    print $0;
	}
    }
}
    
    
