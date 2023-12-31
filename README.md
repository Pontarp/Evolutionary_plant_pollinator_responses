# Evolutionary_plant_pollinator_responses
Evolutionary plant-pollinator responses to anthropogenic land-use change: impacts on ecosystem services

Note: See also Appendix 1 associated with the paper titled "Evolutionary plant-pollinator responses to anthropogenic land-use change: impacts on ecosystem services"


Overview

We analyse an eco-evolutionary model of interacting plants and pollinators to showcase the utility of dynamic models for studying land-use-induced effects on population abundances, selection on traits, and phenotypic adaptation. Trait-based ecological models are well-suited for this purpose, as they provide an abstraction of the key attributes of landscapes and the eco-evolutionary dynamics of interacting organisms within (Loeuille et al. 2013; Georgelin & Loeuille 2016; Pontarp et al. 2019). The landscape is generally modelled as the spatiotemporal distribution of resources, in this particular case two competing plant (R1 and R2) resources on which two competing pollinators (B1 and B2) rely. The two pollinators are assumed to be adapted to the two plants respectively and insect-insect competition is modelled through a trait-based model implementation. We use concepts from niche theory and abstractions of functional traits and trait-matching among populations are assumed to dictate environmental tolerance and ecological interactions. Although generally formulated in the model such matching traits may, for example, be body size which often scales allometrically to many other functional and life history traits in both plants and insects (Niklas 2004; Kalinkat et al. 2015). Alternatively, more specific traits of relevance for environmental preference and ecological interactions including philological traits, phenology, or feeding morphology may be considered (Burkle et al. 2013; Shipley et al. 2017). Regardless of what trait is modelled, adaptive dynamics theory (Geritz et al. 1998; Brännström et al. 2013) lends itself to eco-evolutionary analyses of such trait-based ecological models. Per capita growth (fitness) is trait-dependent and selection on functional traits is thus an emergent property of the landscape properties and the modelled interactions among species. Fitness, in turn, dictates adaptation and the eco-evolutionary dynamics thereof. We utilize these model features in our endeavour to showcase land-use-induced effects on population abundances, selection on traits, and phenotypic adaptation in three essential steps. First, we analyse plant and pollinator population abundances at ecological equilibrium and an evolutionarily stable state (i.e. optimal adaptation in traits). Second, we model land-use change by increasing the carrying capacity for one of the plants (R1) and we analyse ecological time scale changes in abundances. We also analyse change-induced selection pressures on insect traits that dictate trophic mutualistic interactions. Finally, we allow insects to evolve according to the induced selection pressures and we quantify the new ecological equilibrium and evolutionary stable state.  

Ecological model 

We base our model on niche theory and several previous trait-based and eco-evolutionary models (e.g. Christiansen & Loeschcke 1980; Brown & Vincent 1987; Dieckmann & Doebeli 1999). The per-capita growth of two competing plants (R1 and R2) and pollinators (B1 and B2) are formulated as:


(dR_i)/(R_i dt)=g_R  exp⁡(-1/2 ((r_i-r_(m_i ))/σ_r )^2 )-(p_(R_h ) R_i+p_(R_l ) R_j(≠i) )/K_(R_i ) +∑_(j=1,2)▒〖c_(R_j )  exp⁡(-1/2 ((γ_i-β_j)/σ)^2 ) B_j 〗  	(eq. 1)

(dB_i)/(B_i dt)=g_B  exp⁡(-1/2 ((b_i-b_(m_i ))/σ_b )^2 )-(p_(B_h ) B_i+p_(B_l ) B_j(≠i) )/K_(B_i ) +∑_(j=1,2)▒〖c_(B_j )  exp⁡(-1/2 ((β_i-γ_j)/σ)^2 ) R_j 〗 	(eq. 2)


where gR and gB denote the intrinsic growth of plants and pollinators respectively. Similar to the models cited above we adopt a trait-matching approach in modelling environmental tolerance and ecological interactions. More specifically, we model environmental adaptation in plant i through the matching between plant trait ri and the optimal adaptive trait rmi. The plant i is thus fully adapted to the environmental conditions in which it exists when ri = rmi. Adaptation and thus growth rate decreases when ri ≠ rmi and the rate at which such a decrease is modelled as a function of the mismatch mediated through the niche width parameter σr. A large σr models a small decline in growth for a given trait mismatch compared to a small σr which models a larger decline for the same trait mismatch. We model pollinator environmental adaptation in the same way but with different notations for the pollinator environmental trait value bi, optimal pollinator environmental trait bmi, and pollinator environmental niche width parameter σb.

Plants and pollinators are competing within their respective trophic level according to Lotka-Volterra type competition for available resources using interaction coefficients (pRh  and pBh) and the concept of carrying capacity (KRi and KBi). Competition coefficients include intra-specific competition strengths pRh  and pBh which are both scaled to 1 and inter-specific competition pRl  and pBl which is assumed to be > 1 in our model. The interaction coefficients together with the population abundances Ri and Bi thus dictate the total strength of the competition between populations as they compete for some population-specific carrying capacity KRi and KBi. Furthermore, trait-based mutualistic interactions are modelled with a trait-based approach. We assume that the positive effect of mutualistic interactions affects per capita growth rate according to the match between plant and pollinator traits γi and βi respectively. A perfect match between traits renders the largest possible mutualistic effect and similar to the way we model competition we use niche width parameter σ to model the drop-off in effect as a function of trait mismatch. To moderate the number of mutualistic effects that occurs in the system, i.e. the amount of positive effects that flow between trophic levels, conversion coefficients cRj and cBj are used. 

Eco-Evolutionary analysis

We adopt ideas and assumptions from the adaptive dynamics framework (Metz et al. 1992; Geritz et al. 1998; Brännström et al. 2013) to analyse the ecological model described above in an eco-evolutionary context. Such an approach involves four essential components including 1) finding the ecological equilibrium given the ecological model, 2) evaluating invasion fitness of invading mutants in an environment defined by the resident populations at ecological equilibrium, 3) evaluating fitness landscapes and thus potential selection on traits, 4) find the evolutionary stable state (ESS), i.e. when the system reaches optimal adaptation in traits and no more trait changes occur. More specifically, we implement our model described above in R (4.2.0) and we solve for equilibrium population sizes vectors R* and B* by using the ‘solve’ function in R (4.2.0). We use this approach to compute population abundances at ecological equilibrium throughout our analyses and we evaluate the stability of each equilibrium using eigenvalue analysis of the Jacobian matrix of the system. We always compute eigenvalues at the equilibrium point R* and B* and all eigenvalues of the Jacobian matrix had negative real parts throughout our analyses meaning that unstable equilibria are not an issue. Furthermore, focusing on the potential evolution of the pollinator trait β we formulate invasion fitness (per capita growth when rare) of an arbitrary pollinator mutant trait β’ that appears in the environment dictated by the resident trait vectors γ and β, and equilibrium populations size R* and B*: 
	

g_B  exp⁡(-1/2 ((b_i-b_(m_i ))/σ_b )^2 )-(p_(B_h ) B_i+p_(B_l ) B_j(≠i) )/K_(B_i ) +∑_(j=1,2)▒〖c_(B_j )  exp⁡(-1/2 ((β_(i,t)-γ_j)/σ)^2 ) R_j 〗        (i=1,2)		(eq.3)

where β_(i,t) is the invasion species trait values evolving across time. By evaluating the invasion fitness of multiple mutant traits in the same trait dimension as the evolving trait of interest, the fitness landscapes around the focal trait can be quantified. It also follows from the above that the selection gradient can be formulated as the derivative of the invasion fitness function with respect to the focusing traits. When the selection gradient equals zero and the curvature of the selection gradient is concave downwards with a maximum at the resident trait, then ESS is established.

	
Specific model analyses and model parameters   

We start the analyses by setting up the model described above with two plants and two pollinators using the following trait values and parameters: gR  = 1, gB = 1; ri = rmi(i = 1,2), bi = bmi(i = 1,2) (ri and bi are at their optimal),  β1=1.5, β2 =1.5; KR1 = 3, KR2 = 3, KB1 = 1, KB2 = 1; pRh = 1, pRl = 0.5, pBh = 1, pBl = 0.5; γ1 = 1, γ2 = 3; σ = 1; cR1 = 0.1, cR2 = 0.05, cB1 = 0.05, cB2 = 0.1. We then iterate over the essential analysis steps described in the section above. We compute community equilibrium by solving equations 1 and 2. With resident population(s) at equilibrium, we evaluate invasion fitness of β’ values throughout the β trait/ niche space and allow the system to evolve following the canonical equation,
β_(i,t+1)= β_(i,t)+1/2 _^2 s_ 〖f〗_(_i ) (_1,_2 ) B_i          (i=1,2)                     (eq. 4)
where β_(i,t) , 〖f〗_(_i ) (_1,_2 ) and B_i, represent the value of i-th pollinator’s trait, selection gradient, and the population size at the time t correspondingly, ,〖 〗_^2, s_ represents the mutation rate, mutation variance, and mutation steps. Then, we recalculate the ecological equilibrium before progressing iteratively into subsequent evolutionary steps, reiterating the eco-evolutionary procedure described above until the ESS is retrieved, i.e. when all the traits reach their maximum fitness, no mutant trait have positive invasion fitness and no further evolution thus occur. For the model described above and given the parameter set for our analyses, the ESS was found to be at R1=2.18, R2=2.18, B1=0.82, B2=0.82, β1=1.18, β2 =2.82. (see Fig. 2c in the main text). Thereafter we model a disturbance to the system by changing the carrying capacity of plant R1 (KR1) from 3 to 4. We evaluate the new ecological equilibrium and the fitness landscape induced by the change. Population abundances change to R1=1.60, R2=-0.85, B1=0.25, B2=-0.19, and selection is induced on the pollinator traits (see Fig. 2d in the main text). Finally, we analyse the new ESS which renders population abundances R1=4.04, R2=1.04, B1=1.07, B2=0.68, and new traits β1=1.04, β2 =1.18.  (see Fig. 2e in the main text).    

References


Brännström, Å., Johansson, J. & von Festenberg, N. (2013). The hitchhiker’s guide to adaptive dynamics. Games, 4, 304-328.


Brown, J.S. & Vincent, T.L. (1987). A Theory for the Evolutionary Game. Theor Popul Biol, 31, 140-166.


Burkle, L.A., Marlin, J.C. & Knight, T.M. (2013). Plant-Pollinator Interactions over 120 Years: Loss of Species, Co-Occurrence, and Function. Science, 339, 1611-1615.


Christiansen, F.B. & Loeschcke, V. (1980). Evolution and Intraspecific Exploitative Competition .1. One-Locus Theory for Small Additive Gene Effects. Theor Popul Biol, 18, 297-313.


Dieckmann, U. & Doebeli, M. (1999). On the origin of species by sympatric speciation. Nature, 400, 354-357.


Georgelin, E. & Loeuille, N. (2016). Evolutionary response of plant interaction traits to nutrient enrichment modifies the assembly and structure of antagonistic-mutualistic communities. J Ecol, 104, 193-205.

Geritz, S.A.H., Kisdi, E., Meszena, G. & Metz, J.A.J. (1998). Evolutionarily singular strategies and the adaptive growth and branching of the evolutionary tree. Evol Ecol, 12, 35-57.


Kalinkat, G., Jochum, M., Brose, U. & Dell, A.I. (2015). Body size and the behavioral ecology of insects: linking individuals to ecological communities. Curr Opin Insect Sci, 9, 24-30.


Loeuille, N., Barot, S., Georgelin, E., Kylafis, G. & Lavigne, C. (2013). Eco-Evolutionary Dynamics of Agricultural Networks: Implications for Sustainable Management. Advances in Ecological Research, Vol 49: Ecological Networks in an Agricultural World, 49, 339-435.


Metz, J.A.J., Nisbet, R.M. & Geritz, S.A.H. (1992). How should we define fitness for general ecolgical scenarios. Trends Ecol Evol, 7, 198-202.


Niklas, K.J. (2004). Plant allometry: is there a grand unifying theory? Biological Reviews, 79, 871-889.


Pontarp, M., Brännström, A. & Petchey, O.L. (2019). Inferring community assembly processes from macroscopic patterns using dynamic eco-evolutionary models and Approximate Bayesian Computation (ABC). Methods Ecol Evol, 10, 450-460.


Shipley, B., Belluau, M., Kühn, I., Soudzilovskaia, N.A., Bahn, M., Penuelas, J. et al. (2017). Predicting habitat affinities of plant species using commonly measured functional traits. J Veg Sci, 28, 1082-1095.

