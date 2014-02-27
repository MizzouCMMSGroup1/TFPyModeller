# from: http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/ramachandran/calculate/


# import Bio.PDB
# for model in Bio.PDB.PDBParser().get_structure("3U1W_A", "3U1W_A.pdb") :
#     for chain in model :
#         poly = Bio.PDB.Polypeptide.Polypeptide(chain)
#         print "Model %s Chain %s" % (str(model.id), str(chain.id)),
#         print poly.get_phi_psi_list()


import Bio.PDB
for model in Bio.PDB.PDBParser().get_structure("3U1W_A", "3U1W_A.pdb") :
    for chain in model :
        polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides) :
            print "Model %s Chain %s" % (str(model.id), str(chain.id)),
            print "(part %i of %i)" % (poly_index+1, len(polypeptides)),
            print "length %i" % (len(poly)),
            print "from %s%i" % (poly[0].resname, poly[0].id[1]),
            print "to %s%i" % (poly[-1].resname, poly[-1].id[1])
            print poly.get_phi_psi_list()