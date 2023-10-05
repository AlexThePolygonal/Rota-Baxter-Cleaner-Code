from rblib import *


classifier = ActionMnClassifier(2, sg.SR(sg.Integer(0)))
# classifier = ActionMnClassifier(2, sg.SR.var('w_0'))
classifier.classify()

print(f"The number of classes is {len(classifier.QQ_comps)}")

names = ["C1", "C2", "C3", "C4", "A", "B", "A∩B"]
print('Below is the description of the classes')
for i in range(len(classifier.QQ_comps)):
    classifier.QQ_comps[i].name = names[i]
    print(classifier.QQ_comps[i])

print('\n\n\n')

singular_locus = OperatorComponent(classifier.singular_locus, classifier.rota_mat, "A∩B")

e11, e12, e21, e22 = classifier.gen_vars

gubarev_orbits_base = [
    {e21 : e12}, # M1
    {e12 : e21}, # M1 Transpose
    {e21 : e11}, # M2
    {e12 : e11}, # M2 Transpose
    {e21 : e11, e22 : e12}, # M3
    {e12 : e11, e22 : e21}, # M3 Transpose
    {e21 : -e11, e11 : e12}, # M4
    {e12 : -e11, e11 : e21}  # M4 Transpose
]

gubarev_orbits = list(map(classifier.from_gub_repr, gubarev_orbits_base))

for orbit in gubarev_orbits:
    print("The orbit of the point")
    print(orbit)
    print("lies in the component " + (classifier.get_components_of_orbit(classifier.conjugate_by_general_aut(orbit))))
    print("The map from SL_2 is given by")
    print(classifier.conjugate_by_general_aut(orbit))
    print("And the stabilizer by the ideal")
    print((classifier.preimage_of(orbit, classifier.conjugate_by_general_aut(orbit)) + classifier.affine_sl_locus).radical())
    print()

def get_stabiliser(orbit):
    return classifier.preimage_of(orbit, classifier.conjugate_by_general_aut(orbit))



# for orbit in gubarev_orbits:
#     print(classifier.conjugate_by_general_aut(orbit))

# print()

# stabilisers = list(map(lambda stab : (stab + classifier.affine_sl_locus).radical(),map(get_stabiliser, gubarev_orbits)))

# for s in stabilisers:
#     print(s)


