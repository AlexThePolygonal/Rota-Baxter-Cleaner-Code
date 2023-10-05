import sage
import sage.all as sg
from timeit import default_timer as timer
from sage.structure.sage_object import SageObject
from sage.all import singular as singular


_ = singular.LIB("primdec.lib")

def pretty_print_trimmed(ideal, smth):
    gens = sorted(ideal.gens(), key= lambda x : x.degree())
    lins = [gen for gen in gens if gen.degree() == 1]
    substs = {}
    all_vars = set()

    occurs = {}
    for gen in gens:
        for var_ in gen.variables():
            all_vars.add(var_)
            if occurs.get(var_) == None:
                occurs[var_] = 1
            else:
                occurs[var_] += 1
    for lin in lins:
        to_elim = None
        for v in lin.variables():
            if occurs[v] == 1:
                to_elim = v
        if to_elim is not None:
            coeff = lin.monomial_coefficient(to_elim)
            lin_ = lin / coeff - to_elim
            substs[to_elim] = lin_
            gens.remove(lin)
    clean_vars = all_vars - set(substs.keys())
    clean_vars = list(clean_vars)
    vars_long = sg.var(','.join(list(map(chr, range(97, 123)))))
    for v, vv in zip(clean_vars, vars_long):
        substs[v] = vv
    gens = [gen.subs(substs) for gen in gens]
    return smth.subs(substs).subs(substs).subs(substs), gens

def get_entries(mat : sg.Matrix) -> list[sg.RingElement]:
    '''
    Helper func
    '''
    res = []
    for row in mat:
        for entry in row:
            res.append(entry)
    return res

def is_absolutely_prime(sage_QQ_ideal) -> bool:
    '''
    Uses singular to check whether an ideal is prime over the algebraic closure
    '''
    S = singular.ideal(sage_QQ_ideal).absPrimdecGTZ()
    singular.setring(S)
    return singular.execute("size(absolute_primes)") == '1'

def singular_locus(ideal):
    return ideal
    I = sg.macaulay2(ideal, "I")
    J = I.minimalPresentation()
    singlocJ = J.singularLocus().ideal().radical()
    phi = I.dot("cache").dot("minimalPresentationMap")
    return (phi.preimage(singlocJ).radical().sage() + ideal).radical()



class OperatorComponent(SageObject):
    ideal = None
    rota_mat = None
    dimension = -1
    is_absolutely_prime = True
    is_smooth = True
    singular_locus = None
    name = "this class"
    print_minimal_present = True
    
    def __init__(self, ideal, rota_mat, name = "this class"):
        self.ideal = ideal
        self.rota_mat = rota_mat
        self.dimension = ideal.dimension()
        # self.is_absolutely_prime = is_absolutely_prime(ideal)
        self.name = name
        
        self.singular_locus = singular_locus(ideal)
        self.is_smooth = self.singular_locus.dimension() == 0 # the Aff variety is a Proj cone
    
    def contains_orbit(self, orbit):
        for eq in self.ideal.gens():
            if not (eq.subs(orbit.rota_dictionary)) == 0:
                return False
        return True
        

    def _latex_(self):
        '''
        UHHH OHHH bad code
        '''
        print("TODO: Does not work!!!")
        
        I = sg.macaulay2(self.ideal, "I")
        J = I.minimalPresentation().sage()
        Jring = sage.rings.quotient_ring.QuotientRing(J.ring(), J)
        str_vars = [str(gen) for gen in J.ring().gens()]
        base_vars = [str(gen) for gen in self.ideal.ring().gens()]
        str_vars = list(set(base_vars) - set(str_vars))
        str_vars = [gen[0] + '_{' + gen[1:] + "}" for gen in str_vars]
        res_string = self.rota_mat._latex_()
        for var in str_vars:
            res_string = res_string.replace(var, '0')
        res = (res_string + '\n' + Jring._latex_()) 
        return res

    
    def _repr_(self):
        """
        UHHH bad code
        """
        matrix, ideal = pretty_print_trimmed(self.ideal, self.rota_mat)
        res = (f'Operators in {self.name} have the following shape:\n' + matrix.__repr__() + 
            f'\nand their elements must lie in the ideal ' + str(ideal) + ' over the polynomial ring' + 
            f"\nIts dimension as a projective variety is {self.dimension - 1}" + 
            ("\nIt is smooth" if self.is_smooth else f"\nIts singular locus is the {self.singular_locus}" )) + '\n' 
        return res


            
        

def get_components_of_variety(ideal, rota_mat):
    QQ_comps = ideal.minimal_associated_primes()
    sorted(QQ_comps, key=sg.dimension)

    for comp in QQ_comps:
        if not is_absolutely_prime(comp):
            raise NotImplementedError("Non-rational compoments are not supported")
    return list(map(lambda comp : OperatorComponent(comp, rota_mat), QQ_comps))

class RBClassifier(SageObject):
    '''
        Inherit from this class to classify RB-operators on an algebra
        The child class must implement the following:
          `self.gens` --- the generators of the algebra
          `self.to_vector` --- decompose an element using the generators
          `self.to_matrix` --- the inverse to `self.to_vector`
          `self.gen_vars` --- the vars corresponding to gens
    '''

    class OpGubarevRepr():
        gub_dictionary = {}
        rota_dictionary = {}
        rota_mat = None
        ring = sg.SR

        def __init__(self, rbclass, d) -> None:
            # for k, v in d.items():
            #     if not k in rbclass.gen_vars or not v in rbclass.gen_vars:
            #         print(k)
            #         print(v)
            #         print()
            #         print(d)
            #         raise RuntimeError("Error during construction of Rota Operator from Gubarev representation")
            self.gub_dictionary = d
            to_subst = dict([
                (var, 0) for var in rbclass.rota_vars
            ])
            for i in range(rbclass.dim):
                for j in range(rbclass.dim):
                    ei = rbclass.gen_vars[i]
                    ej = rbclass.gen_vars[j]
                    coeff = 0
                    temp = d.get(ei, 0)
                    if temp != 0:
                        coeff = temp.coefficient(ej)
                    to_subst[rbclass.rota_vars[j + i * rbclass.dim]] = coeff
            self.rota_dictionary = to_subst
            self.rota_mat = rbclass.rota_mat.subs(to_subst)
        
        def _latex_(self):
            return self.rota_mat._latex()

        def __repr__(self):
            return self.rota_mat.__repr__()
        
        def _repr_(self):
            return self.rota_mat.__repr__()
        
        def preimage_of(self, pt):
            res = []
            for a, b in zip(self.rota_dictionary.values(), pt.rota_dictionary.values()):
                res.append(a-b)
            return self.ring.ideal(res)
    
    def from_gub_repr(self, d) -> OpGubarevRepr:
        return self.OpGubarevRepr(self, d)


    def __init__(self, weight) -> None:
        '''
        Generate the equations defining the RB operator moduli space
        '''
        self.dim = len(self.gens)
        self.weight = weight
        
        self.rota_vars = sg.SR.var(','.join(f'r{i+1}{j+1}' for i in range(self.dim) for j in range(self.dim)))
        number_of_weight_vars = len(weight.variables())
        self.vars = self.rota_vars + weight.variables()
        self.ring = sage.rings.polynomial.polynomial_ring_constructor.PolynomialRing(sg.QQ, self.vars)
        self.vars = self.ring.gens()
        self.rota_vars = self.vars[:-number_of_weight_vars] if number_of_weight_vars > 0 else self.vars
        
        self.rota_mat = sg.Matrix(sg.SR, self.dim, self.dim)
        for i in range(self.dim):
            for j in range(self.dim):
                self.rota_mat[i,j] = self.rota_vars[i*self.dim + j]

        R = self.rota
        raw_eqns = []
        for x in self.gens:
            for y in self.gens:
                cur_eqns = R(x)*R(y) - (R(R(x)*y + x*R(y)) - weight*R(x*y))
                raw_eqns.append(get_entries(cur_eqns))
        self.raw_eqns = get_entries(raw_eqns)

    def apply_matrix_op(self, op_mat, alg_elem):
        return self.to_matrix(op_mat.T * self.to_vector(alg_elem))
    
    def rota(self, alg_elem):
        '''
        The RB operator
        '''
        return self.apply_matrix_op(self.rota_mat, alg_elem)
         

    def classify(self):
        self.irrelevant_ideal = self.ring.ideal(self.ring.gens())

        self.ideal = self.ring.ideal(self.raw_eqns)
        self.QQ_comps = get_components_of_variety(self.ideal, self.rota_mat) 
       
        self.intersection_locus = self.ring.ideal([1])
        for i in range(len(self.QQ_comps)):
            for j in range(i+1, len(self.QQ_comps)):
                a = self.QQ_comps[i].ideal
                b = self.QQ_comps[j].ideal
                self.intersection_locus = self.intersection_locus.intersection((a+b).radical())
        
        self.singular_locus = self.intersection_locus
        for comp in self.QQ_comps:
            self.singular_locus.intersection(comp.singular_locus)
    

        self.QQ_singloc_comps = get_components_of_variety(self.singular_locus, self.rota_mat)

    def conjugate_by(self, aut, aut_inv, op_mat = None):
        if op_mat is None:
            op_mat = self.rota_mat
        if isinstance(op_mat, self.OpGubarevRepr):
            op_mat = op_mat.rota_mat
        def op(alg_elem):
            return self.apply_matrix_op(op_mat, alg_elem)
        general_mat = sum([gen_var * gen for gen_var, gen in zip(self.gen_vars, self.gens)])
        general_mat_result = aut_inv(op(aut(general_mat)))
        
        def take_var(smth, var):
            subst_dict = dict([
                (gen, 0) for gen in self.gen_vars
            ])
            subst_dict[var] = 1;
            return smth.subs(subst_dict)
        
        gub_dict = dict([
            (ei, (sg.matrix(self.gen_vars) *  self.to_vector(take_var(general_mat_result, ei)))[0][0]) for ei in self.gen_vars
        ])
        return self.OpGubarevRepr(self, gub_dict)
    
    def preimage_of(self, gub_orbit, pt, ring = None):
        if ring is None:
            ring = self.sl_ring
        res = ring.ideal(list([a-b for a, b in [(gub_orbit.rota_dictionary[var], pt.rota_dictionary[var]) for var in self.vars]]))
        return res
    
    def get_components_of_orbit(self, orbit):
        res = []
        for i in range(len(self.QQ_comps)):
            comp = self.QQ_comps[i]
            name = str(i)
            if comp.name != "this class":
                name = comp.name
            if comp.contains_orbit(orbit):
                res.append(name)
        return "∩".join(res)






class MnClassifier(RBClassifier):
    '''
    Classify RB-operators on the matrix algebra Mₙ(ℂ)
    '''
    def __init__(self, dim, weight) -> None:
        self.alg = sg.MatrixSpace(sg.SR, dim)
        self.gens = self.alg.basis().values()
        self.gen_vars = sg.SR.var(','.join([f'e{i+1}{j+1}' for i in range(dim) for j in range(dim)]))
        super().__init__(weight)

    def to_vector(self, alg_elem):
        return sg.Matrix(get_entries(alg_elem)).T
    
    def to_matrix(self, alg_elem_in_basis):
        return sum([self.gens[i] * alg_elem_in_basis[i][0] for i in range(len(self.gens))])
    
    
class ActionMnClassifier(MnClassifier):
    sl_ring = None
    sl_vars = None
    aut = None
    invaut = None
    affine_sl_locus = None
    projective_sl_locus = None


    def __init__(self, dim, weight) -> None:
        super().__init__(dim, weight)
        vars_long = list(map(chr, range(97, 123)))[:dim**2]
        self.sl_vars = sg.SR.var(','.join(vars_long))
        self.sl_ring = sage.rings.polynomial.polynomial_ring_constructor.PolynomialRing(sg.QQ, self.sl_vars)
        self.sl_vars = self.sl_ring.gens()
        phi = sg.Matrix(self.sl_ring, dim, dim)
        for i in range(dim):
            for j in range(dim):
                phi[i,j] = self.sl_vars[i*dim + j]

        phi = phi
        adjphi = phi.adjugate()
        self.projective_sl_locus = self.sl_ring.ideal([phi.det()])
        self.affine_sl_locus = self.sl_ring.ideal([phi.det()-1])
        
        self.aut = lambda x : phi * x * adjphi
        self.invaut = lambda x : adjphi * x * phi


    def conjugate_by_general_aut(self, op_mat):
        res = super().conjugate_by(self.aut, self.invaut, op_mat)
        res.ring = self.sl_ring
        return res