goBambus -C bambus_input.conf -c bambus_input.contig -m bambus_input.mates -x bambus_input.xml -o bambus_out
dot -Tpdf -o bambus_out.pdf bambus_out.dot
untangle -e bambus_out.evidence.xml -s bambus_out.out.xml -o scaffold.untangle.xml
printScaff -e bambus_out.evidence.xml -s scaffold.untangle.xml -l bambus_out.lib -dot -detail -oo -sum -f nucmer_input.contig -nomerge -o scaffold
dot -Tpdf -o scaffold.pdf scaffold.dot

