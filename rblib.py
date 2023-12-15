import sage
import sage.all as sg
from timeit import default_timer as timer
from sage.structure.sage_object import SageObject
from sage.all import singular as singular
import copy
import warnings
import random
import itertools
import typing

random.seed(441)
sg.set_random_seed(442)
sg.macaulay2.set_seed(443)

_ = singular.LIB("primdec.lib")



def batched_it(iterable, n):
    "Batch data into iterators of length n. The last batch may be shorter."
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while True:
        chunk_it = (itertools.islice(it, n))
        try:
            first_el = next(chunk_it)
        except StopIteration:
            return
        yield itertools.chain((first_el,), chunk_it)

def print_eqn_system(
        ideal_gens : list[sage.structure.element.RingElement], 
        eq_per_line : int = 3
    ):
    """
    Print the equations in the align* environment
    """
    gens = iter(map(sg.latex, list(ideal_gens)))
    text = "\\begin{align*}\n" + " = 0 \\\\\n".join([
        " = 0 \\quad\n".join(chunk) for chunk in batched_it(gens, eq_per_line)
    ]) +  " = 0 \n\\end{align*}"
    return text

def get_ideal_simplification(
        ideal : sage.rings.ideal.Ideal
    ) -> tuple[
        dict[sage.structure.element.RingElement, sage.structure.element.RingElement], 
        list[sage.structure.element.RingElement]
        ]:
    """
    Simplify the given ideal, eliminating the linearly dependent variables
    """
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
            substs[to_elim] = -lin_
            gens.remove(lin)
    clean_vars = all_vars - set(substs.keys())
    clean_vars = list(clean_vars)
    vars_long = sg.var(','.join(list(map(chr, range(97, 123)))))
    for v, vv in zip(clean_vars, vars_long):
        substs[v] = vv
    gens = [gen.subs(substs) for gen in gens]
    return substs, gens
 

def pretty_print_trimmed(ideal : sage.rings.ideal.Ideal, smth):
    substs, gens = get_ideal_simplification(ideal)
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

def is_absolutely_prime(sage_QQ_ideal : sage.rings.ideal.Ideal) -> bool:
    '''
    Uses singular to check whether an ideal is prime over the algebraic closure
    '''
    S = singular.ideal(sage_QQ_ideal).absPrimdecGTZ()
    singular.setring(S)
    return singular.execute("size(absolute_primes)") == '1'

def singular_locus(ideal : sage.rings.ideal.Ideal) -> sage.rings.ideal.Ideal:
    # return ideal
    I = sg.macaulay2(ideal, "I")
    J = I.minimalPresentation()
    singlocJ = J.singularLocus().ideal().radical()
    phi = I.dot("cache").dot("minimalPresentationMap")
    return (phi.preimage(singlocJ).radical().sage() + ideal).radical()



class OperatorComponent(SageObject):
    ideal : sage.rings.ideal.Ideal = None
    rota_mat : sage.matrix.constructor.Matrix = None
    dimension : int = -1
    is_absolutely_prime : bool = True
    is_smooth : bool = True
    singular_locus : typing.Optional[sage.rings.ideal.Ideal] = None
    name : str = "this class"
    print_minimal_present : bool = True
    
    def __init__(self, ideal : sage.rings.ideal.Ideal, rota_mat : sage.matrix.constructor.Matrix, name : str = "this class") -> None:
        self.ideal = ideal
        self.rota_mat = rota_mat
        self.dimension = ideal.dimension()
        # self.is_absolutely_prime = is_absolutely_prime(ideal)
        self.name = name
        
        self.singular_locus = singular_locus(ideal)
        self.is_smooth = self.singular_locus.dimension() == 0 # the Aff variety is a Proj cone
    
    def contains_orbit(self, orbit : 'RBClassifier.OpRepr') -> bool:
        for eq in self.ideal.gens():
            if not (eq.subs(orbit.rota_dictionary)) == 0:
                return False
        return True
       

    def _latex_(self) -> str:
        '''
        Bad code to print as latex
        '''
        # print("TODO: Does not work!!!") 
        matrix, ideal = pretty_print_trimmed(self.ideal, self.rota_mat)
        res = (f'Operators in {self.name} have the following shape:\n$$' + matrix._latex_() + "$$"
            f'\nand their elements must lie in the ideal ' + print_eqn_system(ideal, 3) + '' + 
            f"\nIts dimension as a projective variety is {self.dimension - 1}." + 
            ("\nIt is smooth" if self.is_smooth else f"\nIts singular locus is the ideal ${self.singular_locus._latex_()}$" )) + '.\n' 
        return res
 

    
    def _repr_(self) -> str:
        """
        bad code to humanly-legible print
        """
        matrix, ideal = pretty_print_trimmed(self.ideal, self.rota_mat)
        res = (f'Operators in {self.name} have the following shape:\n' + matrix.__repr__() + 
            f'\nand their elements must lie in the ideal ' + str(ideal) + ' over the polynomial ring' + 
            f"\nIts dimension as a projective variety is {self.dimension - 1}" + 
            ("\nIt is smooth" if self.is_smooth else f"\nIts singular locus is the {self.singular_locus}" )) + '\n' 
        return res

        

def get_components_of_variety(ideal : sage.rings.ideal.Ideal, rota_mat : sage.matrix.constructor.Matrix) -> list[OperatorComponent]:
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
    gens : sage.structure.element.RingElement
    

    class OpRepr(SageObject):
        """
        Representative of an operator class
        """
        rota_dictionary : dict[sage.structure.element.RingElement, sage.structure.element.RingElement] = {}
        rota_mat : sage.matrix.constructor.Matrix = None
        ring : sage.rings.ring.Ring = sg.SR
        name : str= ""

    class OpGubarevRepr(OpRepr):
        """
        Representative of an operator class in the notation used by Gubarev
        """
        gub_dictionary : dict[sage.structure.element.RingElement, sage.structure.element.RingElement] = {}

        def __init__(self, rbclass, d) -> None:
            # for k, v in d.items():
            #     if not k in rbclass.gen_vars or not v in rbclass.gen_vars:
            #         print(k)
            #         print(v)
            #         print()
            #         print(d)
            #         raise RuntimeError("Error during construction of Rota Operator from Gubarev representation")
            super().__init__()
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
        
        def _latex_(self) -> str:
            return self.rota_mat._latex_()

        def __repr__(self) -> str:
            return self.rota_mat.__repr__()
        
        def _repr_(self) -> str:
            return self.rota_mat.__repr__()
        
        def preimage_of(self, pt):
            res = []
            for a, b in zip(self.rota_dictionary.values(), pt.rota_dictionary.values()):
                res.append(a-b)
            return self.ring.ideal(res)
    
    def from_gub_repr(self, d : dict, name = "") -> OpGubarevRepr:
        """
        Create an OpRepr from orbit written in Gubarev Notation
        """
        res = self.OpGubarevRepr(self, d)
        res.name = ""
        return res
    
    def from_mat(self, mat : sage.matrix.constructor.Matrix, name = "") -> OpRepr:
        res = self.OpRepr()
        res.name = name
        res.rota_mat = mat
        res.rota_dictionary = {}
        for i in range(self.dim):
            for j in range(self.dim):
                res.rota_dictionary[self.rota_vars[i*self.dim + j]] = res.rota_mat[i,j]
        return res
    
    def find_orbit_parametrization(self, orbit : OpRepr, component : OperatorComponent) -> dict[sage.structure.element.RingElement, sage.structure.element.RingElement]:
        warnings.warn("`RBClassifier.find_orbit_parametrization` is very crude. Verify that all variables are represented!")
        substs, gens = get_ideal_simplification(component.ideal)
        orbit_dict : dict = orbit.rota_dictionary
        newvars = set([var for gen in gens for var in gen.variables()])
        temp = dict([(var, None) for var in self.rota_vars])

        for rv in self.rota_vars:
            u = substs[rv]
            if u in newvars:
                temp[rv] = u

        res = dict()
        for r_var, sl_expr in orbit_dict.items():
            if temp[r_var] is not None:
                res[temp[r_var]] = sl_expr

        # raise KeyboardInterrupt()
        return res
    
    def print_orbit(self, orbit : OpRepr, comp_dict : dict[str, OperatorComponent]):
        """
        Prints the orbit laconically if it can find its component in comp_dict
        """
        name = self.get_components_of_orbit(orbit)
        comp = comp_dict[name]
        res = self.find_orbit_parametrization(orbit, comp)
        items = sorted(res.items(), key = lambda x : x[0]._latex_())
        text = '\\begin{align*}\n' + "\\\\\n".join(
            [var._latex_() + " = " + expr.full_simplify()._latex_() for var, expr in items]
        ) + " = 0\n\\end{align*}"
        return text

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

    def apply_matrix_op(self, op_mat : sage.matrix.constructor.Matrix, alg_elem : sage.structure.element.RingElement):
        return self.to_matrix(op_mat.T * self.to_vector(alg_elem))
    
    def rota(self, alg_elem : sage.structure.element.RingElement) -> sage.structure.element.RingElement:
        '''
        The RB operator
        '''
        return self.apply_matrix_op(self.rota_mat, alg_elem)
         

    def classify(self) -> None:
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


    def rewrite_as_rota(self, func):
        v = sum([self.gen_vars[i] * self.gens[i] for i in range(len(self.gens))])
        def take_ith(i):
            substs = dict([(gv, 0) for gv in self.gen_vars])
            substs[self.gen_vars[i]] = 1
            return substs
        matv = self.to_vector(func(v))
        res = self.rota_mat
        def take(i):
            res = dict([(g, 0) for g in self.gen_vars])
            res[self.gen_vars[i]] = 1
            return res
        
        for i in range(self.dim):
            for j in range(self.dim):
                res[i,j] = matv[i].subs(take(j))[0]
        return res




class MnClassifier(RBClassifier):
    alg = sg.MatrixSpace(sg.SR, 1)
    gens = alg.basis().values()
    gen_vars = sg.SR.var("erm,sorry")
    det = 1
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
        # vars_long = list(map(chr, range(97, 123)))[:dim**2]
        vars_long = ['x', 'y', 'z', 't']#, 'w', 'u', 'r', 's', 'p', 'q']
        self.sl_vars = sg.SR.var(','.join(vars_long))
        self.sl_ring = sage.rings.polynomial.polynomial_ring_constructor.PolynomialRing(sg.QQ, self.sl_vars)
        self.sl_vars = self.sl_ring.gens()
        phi = sg.Matrix(self.sl_ring, dim, dim)
        for i in range(dim):
            for j in range(dim):
                phi[i,j] = self.sl_vars[i*dim + j]

        phi = phi
        adjphi = phi.adjugate()
        self.det = sg.det(phi)
        self.projective_sl_locus = self.sl_ring.ideal([phi.det()])
        self.affine_sl_locus = self.sl_ring.ideal([phi.det()-1])
        
        self.aut = lambda x : phi * x * adjphi
        self.invaut = lambda x : adjphi * x * phi

        self.phi = phi
        self.adjphi = adjphi
    

    def conjugate_by_general_aut(self, op_mat):
        if isinstance(op_mat, self.OpGubarevRepr):
            op_mat = op_mat.rota_mat
        def func(x):
            return self.invaut(self.to_matrix(
                    op_mat * self.to_vector(
                        self.aut(x)
                    )
                ))
        mat = self.rewrite_as_rota(func)
        res = self.from_mat(mat)
        res.ring = self.sl_ring
        return res