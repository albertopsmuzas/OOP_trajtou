#!/usr/bin/awk -f
# ******************************************************************************************
#
# USAGE:
# trajtouCountQuantizedTrajs.awk n=? m=? V=? J=? mJ=? filename
#
# ******************************************************************************************
#
# - filename: should be a file with structure like OUTANA6Dmappingpeaks.out
# - n=?  : ? should be an integer or keyword "all". Counts over diffraction states.
# - m=?  : ? should be an integer or keyword "all". Counts over diffraction states.
# - V=?  : ? should be an integer or keyword "all". Counts over vibrational states.
# - J=?  : ? should be an integer or keyword "all". Counts over rotational states. 
# - mJ=? : ? should be an integer or keyword "all". Counts over polarized rotational states.
#
# WARNING: IMPLEMENTED CASES:
# - Only one parameter (It does not matter which one) is initialized with "all" value
# - All parameters but n and m are initialized with "all" value
# - All parameters are initialized with "all" value
# ******************************************************************************************
BEGIN{
	colN=4;
	colM=5;
	colV=6;
	colJ=7;
	colmJ=8;
	count=0;
	condErr=0;
	n="undefined";
	m="undefined";
	V="undefined";
	J="undefined";
	mJ="undefined";
}
!/\#/{
	if ( n == "undefined" || m == "undefined" || V == "undefined" || J == "undefined" || mJ == "undefined" )
	{
		print "At least one parameter was not initialized and still has the value \"undefined\".";
		exit 1;
	}
	else if ( n != "all" && m != "all" && V != "all" && J != "all" && mJ != "all" )
	{
		if ($colN==n && $colM==m && $colV==V && $colJ==J && $colmJ==mJ) count=count+1;
	}
	else if ( n == "all" && m != "all" && V != "all" && J != "all" && mJ != "all" )
	{
		if ($colM=m && $colV==V && $colJ==J && $colmJ==mJ) count=count+1;
	}
	else if ( n != "all" && m == "all" && V != "all" && J != "all" && mJ != "all" )
	{
		if ($colN==n && $colV==V && $colJ==J && $colmJ==mJ) count=count+1;
	}
	else if ( n != "all" && m != "all" && V == "all" && J != "all" && mJ != "all" )
	{
		if ($colN==n && $colM==m && $colJ==J && $colmJ==mJ) count=count+1;
	}
	else if ( n != "all" && m != "all" && V != "all" && J == "all" && mJ != "all" )
	{
		if ($colN==n && $colM==m && $colV==V && $colmJ==mJ) count=count+1;
	}
	else if ( n != "all" && m != "all" && V != "all" && J != "all" && mJ == "all" )
	{
		if ($colN==n && $colM==m && $colV==V && $colJ==J) count=count+1;
	}
	else if ( n != "all" && m != "all" && V == "all" && J == "all" && mJ == "all" )
	{
		if ($colN==n && $colM==m) count=count+1;
	}
	else if ( n == "all" && m == "all" && V != "all" && J != "all" && mJ != "all" )
	{
		if ($colV==V && $colJ==J && $colmJ==mJ) count=count+1;
	}
	else if ( n == "all" && m == "all" && V == "all" && J == "all" && mJ == "all" )
	{
		count=count+1;
	}
	else
	{
		print "Specification not implemented. Check script";
		exit 1;
	}
}
END{
		print count;
}
