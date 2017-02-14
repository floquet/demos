#!/bin/csh
foreach i (*.for)
set bname=`basename $i ".for" `
mv $bname.for $bname.f
end

