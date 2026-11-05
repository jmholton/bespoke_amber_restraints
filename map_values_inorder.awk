#! /bin/awk -f
#
#  put phenix.map_value_at_point output back in the order of the input atoms     -James Holton 11-5-26
#
#  phenix reads the PDB into a hierarchy and prints the points in ITS order: within an amino
#  acid that is N CA C O CB ..., not Amber's N CA CB ... C O.  So never count "Map value:" lines.
#  Usage:  cat probe.pdb phenix.log | map_values_inorder.awk
#  prints one "Map value: <rho>" line per ATOM/HETATM of probe.pdb, in its order, matched by the
#  printed coordinates (0 and a warning on stderr if phenix gave none for an atom)
#
/^ATOM|^HETAT/{key[++n]=sprintf("(%.3f,%.3f,%.3f)",substr($0,31,8),substr($0,39,8),substr($0,47,8));next}
# "Input point: (x,y,z) Fractionalized: (...) Map value: v"
/Map value:/ && $3~/^\(/{rho[$3]=$NF;next}
# older phenix: "..." Point: x, y, z ... Map value: v
/Map value:/ && /^"/{point=substr($0,index($0," Point: "));split(point,w);rho["(" w[2] w[3] w[4] ")"]=$NF;next}
END{for(i=1;i<=n;++i){
      if(key[i] in rho) print "Map value:",rho[key[i]];
      else {print "Map value: 0";++miss}}
    if(miss) printf("WARNING: map_values_inorder: no map value for %d of %d atoms\n",miss,n) > "/dev/stderr"}
