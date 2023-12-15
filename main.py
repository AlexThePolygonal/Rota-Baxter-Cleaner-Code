from rblib import *


classifier = ActionMnClassifier(2, sg.SR(sg.Integer(0)))
# classifier = ActionMnClassifier(2, sg.SR.var('w_0'))
classifier.classify()

print(f"The number of classes is ${len(classifier.QQ_comps)}$.")

names = ["C1", "C2", "C3", "C4", "A", "B", "A∩B"]
print('The projective variety of Rota-Baxter operators consists of 6 components: four curves $C_1, C_2, C_3, C_4$ and two surfaces $A$ and $B$.')
print("\\begin{itemize}")
for i in range(len(classifier.QQ_comps)):
    print("\item")
    classifier.QQ_comps[i].name = names[i]
    print((classifier.QQ_comps[i])._latex_())
    # print("----------")
    # print(classifier.QQ_comps[i].rota_mat)
    # print(classifier.QQ_comps[i].ideal)
    # print()
print("\\end{itemize}")
print('\n\n\n')

singular_locus = OperatorComponent(classifier.singular_locus, classifier.rota_mat, "A∩B")

print(singular_locus._latex_())

total_comps = dict([(comp.name, comp) for comp in (classifier.QQ_comps + [singular_locus])])

e11, e12, e21, e22 = classifier.gen_vars

gubarev_orbits_base = [
    ({e21 : e12}, "M_1"),
    ({e12 : e21}, "M_1^T)"),
    ({e21 : e11}, "M_2"),
    ({e12 : e11}, "M_2^T"),
    ({e21 : e11, e22 : e12}, "M_3"),
    ({e12 : e11, e22 : e21}, "M_3^T"),
    ({e21 : -e11, e11 : e12}, "M_4"),
    ({e12 : -e11, e11 : e21}, "M_4^T")
]

gubarev_orbits = [classifier.from_gub_repr(d, name) for d, name in  gubarev_orbits_base]

for orbit in gubarev_orbits:

    phi = classifier.conjugate_by_general_aut(orbit)
    
    print("The orbit of the point")
    print("$$" + orbit._latex_() + "$$")
    print("lies in the component " + (classifier.get_components_of_orbit(phi)))
    print("The map from SL_2 is given by")
    # print(orbit.rota_mat.subs(phi.rota_dictionary))
    print(classifier.print_orbit(phi, total_comps))
    # print("$$" + phi.rota_mat._latex_() + "$$")
    print("And the stabilizer by the ideal")
    print("$" + (classifier.preimage_of(orbit, phi) + classifier.affine_sl_locus).radical()._latex_() + "$")
    
    # for eqn in classifier.raw_eqns:
    #     if (eqn.subs(phi.rota_dictionary) != 0):
    #         print(eqn.subs(phi.rota_dictionary).full_simplify())


    print()


def get_stabiliser(orbit):
    return classifier.preimage_of(orbit, classifier.conjugate_by_general_aut(orbit))



a,b,c = sg.var("a,b,c")
tang_classes_base = [
    ("R_1", sg.matrix([
        [ 0,  0,  1,  a],
        [-1, -a,  0,  0],
        [-0,  0,  a, a*a],
        [-a,-a*a,  0,  0]
    ])),
    ("R_2", sg.matrix([
        [ a,-a*a,-0,  0],
        [ 1, -a,  0,  0],
        [ 0,  0,  a,  -a*a],
        [ 0,  0,  1,  -a]
    ])),
    ("R_3", sg.matrix([
        [ a, -a*a, -0,  0],
        [ 1, -a,  0,  0],
        [-a*a, a**3, -0,  0],
        [ -a,  a*a,  0,  0]
    ])), 
    ("R_4", sg.matrix([
        [ a,  0,  1,  0],
        [ 0,  a,  0,  1],
        [-a*a,0, -a,  0],
        [ 0,-a*a, 0, -a]
    ])),
    ("R_5", sg.matrix([
        [ 0, a*a,-0,  a],
        [ 0,  a,  0,  1],
        [-a*a,0, -a,  0],
        [-a,  0, -1,  0]
    ])),
    ("R_6", sg.matrix([
        [ 0,  1,  0,  a],
        [ 0,  a,  0,a*a],
        [-0,-1/a, 0, -1],
        [ 0, -1,  0, -a]
    ])),
    ("R_7", sg.matrix([
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  1],
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  0]
    ])),
    ("R_8", sg.matrix([
        [ 0,  1,  0,  0],
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  1],
        [ 0,  0,  0,  0]
    ])),
    ("R_9", sg.matrix([
        [ 0,  1,  0,  0],
        [ 0,  0,  0,  0],
        [-1,  0,  0,  0],
        [ 0,  0,  0,  0]
    ])),
    ("R_10", sg.matrix([
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  0],
        [ a,  1,  0,  0],
        [ 0,  0,  0,  0]
    ])),
    ("R_11", sg.matrix([
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  0],
        [ 0,  1,  0,  a],
        [ 0,  0,  0,  0]
    ])),
    ("R_12", sg.matrix([
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  1],
        [ 0,  0,  0,  0]
    ])),
    ("R_13", sg.matrix([
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  1],
        [ 0, -1,  0,  0]
    ])),
    ("R_14", sg.matrix([
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  0],
        [ a,  0,  0,  0],
        [ 0,  0,  0,  0]
    ])),
    ("R_15", sg.matrix([
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  0],
        [ 1,  0,  0,  0],
        [ 0,  1,  0,  0]
    ])),
    ("R_16", sg.matrix([
        [ 0,  0,  a,  1],
        [ 0,  0,-a*a,-a],
        [ 0,  0,  1, 1/a],
        [ 0,  0, -a, -1]
    ])),
    ("R_17", sg.matrix([
        [ 1,  0,  a,  0],
        [ a,  0, a*a, 0],
        [-1/a,0, -1,  0],
        [-1,  0, -a,  0]
    ])),
    ("R_18", sg.matrix([
        [ 0,  0,  0,  0],
        [ 0,  0,  1,  a],
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  0]
    ])),
    ("R_19", sg.matrix([
        [ 0,  0,  0,  0],
        [ a,  0,  1,  0],
        [ 0,  0,  0,  0],
        [ 0,  0,  0,  0]
    ])),
    ("R_20", sg.matrix([
        [ b*b,  a,   b,   a/b],

        [  b,  a/b,  1,  a/(b*b)],

        [-b**3,-a*b,-b*b,  -a],

        [-b*b,  -a,   -b, -a/b]
    ])),
    ("R_21", sg.matrix([
        [  -a/b,  a,   b,  -b*b],

        [-a/(b*b),a/b,  1,  -b],

        [    a,  -a*b,-b*b,b**3],

        [   a/b, -a,   -b, b*b]
    ]))
]

tang_classes = [classifier.from_mat(mat, name) for name, mat in tang_classes_base]


print("Classes from Tang et al(2014) correspond to our classes as follows:")
for tang_class in tang_classes:
    print(tang_class.name + ": " + classifier.get_components_of_orbit(tang_class))

print()

# for orbit in gubarev_orbits:
#     print(classifier.conjugate_by_general_aut(orbit))

# print()

# stabilisers = list(map(lambda stab : (stab + classifier.affine_sl_locus).radical(),map(get_stabiliser, gubarev_orbits)))

# for s in stabilisers:
#     print(s)


