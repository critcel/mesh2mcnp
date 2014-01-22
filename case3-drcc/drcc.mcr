#!MC 900
$!VarSet |num_z|=  2
$!VarSet |cm_per_zlev|=  25
$!FIELDLAYERS SHOWMESH = NO
$!FIELDLAYERS SHOWBOUNDARY = NO
$!FRAMEMODE = THREED
$!GLOBALCONTOUR VAR = 4
 $!VarSet |cm_tot|=(|cm_per_zlev|*|num_z|)
 $!FIELD [1-|cm_tot|] SCATTER{COLOR = MULTI}
 $!FIELD [1-|cm_tot|] SCATTER{ISFILLED = YES}
 $!FIELD [1-|cm_tot|] SCATTER{FILLCOLOR = MULTI}
 $!FIELD [1-|cm_tot|] SCATTER{FRAMESIZE = 1}
 $!LOOP |num_z|
 $!VARSET |a| = ((|loop|-1)*|cm_per_zlev|+1)
 $!VARSET |b| = (|loop|*|cm_per_zlev|)
 $!FIELD [|a|-|b|]  GROUP = |loop|
 $!ENDLOOP
$!FIELDLAYERS SHOWSCATTER = YES
