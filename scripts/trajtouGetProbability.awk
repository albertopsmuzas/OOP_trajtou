#!/usr/bin/awk -f
# ******************************************************************************************
#
# USAGE:
# trajtouCountQuantizedTrajs.awk n=? m=? V=? J=? mJ=? filename
#
# ******************************************************************************************
#
# - filename: should be a file with structure like OUTANA6DseenPeaks.out
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
	colDiffOrder=1;
	colN=2;
	colM=3;
	colV=4;
	colJ=5;
	colmJ=6;
	colProb=12;
	prob=0;
	condErr=0;
	badness=0;
	diffOrder="unset";
	n="undefined";
	m="undefined";
	V="undefined";
	J="undefined";
	mJ="undefined";
}
#!/\#/{
{	if ( $1 != "#" )
	{
		if ( n == "undefined" || m == "undefined" || V == "undefined" || J == "undefined" || mJ == "undefined" )
		{
			print "At least one parameter was not initialized and still has the value \"undefined\".";
			badness=1;
			exit 1;
		}
		else if ( diffOrder =="unset" && n != "all" && m != "all" && V != "all" && J != "all" && mJ != "all" )
		{
			if ($colN==n && $colM==m && $colV==V && $colJ==J && $colmJ==mJ) prob=prob+$colProb;
		}
		else if ( diffOrder =="unset" && n == "all" && m != "all" && V != "all" && J != "all" && mJ != "all" )
		{
			if ($colM=m && $colV==V && $colJ==J && $colmJ==mJ) prob=prob+$colProb;
		}
		else if ( diffOrder =="unset" && n != "all" && m == "all" && V != "all" && J != "all" && mJ != "all" )
		{
			if ($colN==n && $colV==V && $colJ==J && $colmJ==mJ) prob=prob+$colProb;
		}
		else if ( diffOrder =="unset" && n != "all" && m != "all" && V == "all" && J != "all" && mJ != "all" )
		{
			if ($colN==n && $colM==m && $colJ==J && $colmJ==mJ) prob=prob+$colProb;
		}
		else if ( diffOrder =="unset" && n != "all" && m != "all" && V != "all" && J == "all" && mJ != "all" )
		{
			if ($colN==n && $colM==m && $colV==V && $colmJ==mJ) prob=prob+$colProb;
		}
		else if ( diffOrder =="unset" && n != "all" && m != "all" && V != "all" && J != "all" && mJ == "all" )
		{
			if ($colN==n && $colM==m && $colV==V && $colJ==J) prob=prob+$colProb;
		}
		else if ( diffOrder =="unset" && n != "all" && m != "all" && V == "all" && J == "all" && mJ == "all" )
		{
			if ($colN==n && $colM==m) prob=prob+$colProb;
		}
		else if ( diffOrder =="unset" && n == "all" && m == "all" && V != "all" && J != "all" && mJ != "all" )
		{
			if ($colV==V && $colJ==J && $colmJ==mJ) prob=prob+$colProb;
		}
		else if ( diffOrder =="unset" && n == "all" && m == "all" && V == "all" && J != "all" && mJ == "all" )
		{
			if ($colJ==J) prob=prob+$colProb;
		}
		else if ( diffOrder =="unset" && n == "all" && m == "all" && V == "all" && J == "all" && mJ == "all" )
		{
			prob=prob+$colProb;
		}
		else if ( diffOrder != "unset" )
		{
			if($colDiffOrder==diffOrder) prob=prob+$colProb;
		}
		else
		{
			print "Specification not implemented. Check script";
			badness=1;
			exit 1;
		}
	}
	else
	{
		if( $2 == "Total" && $4 == "scattered" ) ntrajs=$11;
	}
}
END{
		if(badness==0) print prob,"+/-",sqrt(prob*(1-prob)/(ntrajs-1));
}
