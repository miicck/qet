import os
base_dir = os.path.dirname(os.path.abspath(__file__))

titles = [
"atomic number", 
"name", 
"symbol", 
"mass number", 
"neutrons", 
"radioactive",
"covalent radius", # https://doi.org/10.1039%2Fb801115j
]

data = [
[1  ,"Hydrogen     ","H ",1.007  ,0  ,0, 0.31],
[2  ,"Helium       ","He",4.002  ,2  ,0, 0.28],
[3  ,"Lithium      ","Li",6.941  ,4  ,0, 1.28],
[4  ,"Beryllium    ","Be",9.012  ,5  ,0, 0.96],
[5  ,"Boron        ","B ",10.811 ,6  ,0, 0.84],
[6  ,"Carbon       ","C ",12.011 ,6  ,0, 0.76],
[7  ,"Nitrogen     ","N ",14.007 ,7  ,0, 0.71],
[8  ,"Oxygen       ","O ",15.999 ,8  ,0, 0.66],
[9  ,"Fluorine     ","F ",18.998 ,10 ,0, 0.57],
[10 ,"Neon         ","Ne",20.18  ,10 ,0, 0.58],
[11 ,"Sodium       ","Na",22.99  ,12 ,0, 1.66],
[12 ,"Magnesium    ","Mg",24.305 ,12 ,0, 1.41],
[13 ,"Aluminum     ","Al",26.982 ,14 ,0, 1.21],
[14 ,"Silicon      ","Si",28.086 ,14 ,0, 1.11],
[15 ,"Phosphorus   ","P ",30.974 ,16 ,0, 1.07],
[16 ,"Sulfur       ","S ",32.065 ,16 ,0, 1.05],
[17 ,"Chlorine     ","Cl",35.453 ,18 ,0, 1.02],
[18 ,"Argon        ","Ar",39.948 ,22 ,0, 1.06],
[19 ,"Potassium    ","K ",39.098 ,20 ,0, 2.03],
[20 ,"Calcium      ","Ca",40.078 ,20 ,0, 1.76],
[21 ,"Scandium     ","Sc",44.956 ,24 ,0, 1.70],
[22 ,"Titanium     ","Ti",47.867 ,26 ,0, 1.60],
[23 ,"Vanadium     ","V ",50.942 ,28 ,0, 1.53],
[24 ,"Chromium     ","Cr",51.996 ,28 ,0, 1.39],
[25 ,"Manganese    ","Mn",54.938 ,30 ,0, 1.39],
[26 ,"Iron         ","Fe",55.845 ,30 ,0, 1.32],
[27 ,"Cobalt       ","Co",58.933 ,32 ,0, 1.26],
[28 ,"Nickel       ","Ni",58.693 ,31 ,0, 1.24],
[29 ,"Copper       ","Cu",63.546 ,35 ,0, 1.32],
[30 ,"Zinc         ","Zn",65.38  ,35 ,0, 1.22],
[31 ,"Gallium      ","Ga",69.723 ,39 ,0, 1.22],
[32 ,"Germanium    ","Ge",72.64  ,41 ,0, 1.20],
[33 ,"Arsenic      ","As",74.922 ,42 ,0, 1.19],
[34 ,"Selenium     ","Se",78.96  ,45 ,0, 1.20],
[35 ,"Bromine      ","Br",79.904 ,45 ,0, 1.20],
[36 ,"Krypton      ","Kr",83.798 ,48 ,0, 1.16],
[37 ,"Rubidium     ","Rb",85.468 ,48 ,0, 2.20],
[38 ,"Strontium    ","Sr",87.62  ,50 ,0, 1.95],
[39 ,"Yttrium      ","Y ",88.906 ,50 ,0, 1.90],
[40 ,"Zirconium    ","Zr",91.224 ,51 ,0, 1.64],
[41 ,"Niobium      ","Nb",92.906 ,52 ,0, 1.64],
[42 ,"Molybdenum   ","Mo",95.96  ,54 ,0, 1.54],
[43 ,"Technetium   ","Tc",98     ,55 ,1, 1.47],
[44 ,"Ruthenium    ","Ru",101.07 ,57 ,0, 1.46],
[45 ,"Rhodium      ","Rh",102.906,58 ,0, 1.42],
[46 ,"Palladium    ","Pd",106.42 ,60 ,0, 1.39],
[47 ,"Silver       ","Ag",107.868,61 ,0, 1.45],
[48 ,"Cadmium      ","Cd",112.411,64 ,0, 1.44],
[49 ,"Indium       ","In",114.818,66 ,0, 1.42],
[50 ,"Tin          ","Sn",118.71 ,69 ,0, 1.39],
[51 ,"Antimony     ","Sb",121.76 ,71 ,0, 1.39],
[52 ,"Tellurium    ","Te",127.6  ,76 ,0, 1.38],
[53 ,"Iodine       ","I ",126.904,74 ,0, 1.39],
[54 ,"Xenon        ","Xe",131.293,77 ,0, 1.40],
[55 ,"Cesium       ","Cs",132.905,78 ,0, 2.44],
[56 ,"Barium       ","Ba",137.327,81 ,0, 2.15],
[57 ,"Lanthanum    ","La",138.905,82 ,0, 2.07],
[58 ,"Cerium       ","Ce",140.116,82 ,0, 2.04],
[59 ,"Praseodymium ","Pr",140.908,82 ,0, 2.03],
[60 ,"Neodymium    ","Nd",144.242,84 ,0, 2.01],
[61 ,"Promethium   ","Pm",145    ,84 ,1, 1.99],
[62 ,"Samarium     ","Sm",150.36 ,88 ,0, 1.98],
[63 ,"Europium     ","Eu",151.964,89 ,0, 1.98],
[64 ,"Gadolinium   ","Gd",157.25 ,93 ,0, 1.96],
[65 ,"Terbium      ","Tb",158.925,94 ,0, 1.94],
[66 ,"Dysprosium   ","Dy",162.5  ,97 ,0, 1.92],
[67 ,"Holmium      ","Ho",164.93 ,98 ,0, 1.92],
[68 ,"Erbium       ","Er",167.259,99 ,0, 1.89],
[69 ,"Thulium      ","Tm",168.934,100,0, 1.90],
[70 ,"Ytterbium    ","Yb",173.054,103,0, 1.87],
[71 ,"Lutetium     ","Lu",174.967,104,0, 1.75],
[72 ,"Hafnium      ","Hf",178.49 ,106,0, 1.87],
[73 ,"Tantalum     ","Ta",180.948,108,0, 1.70],
[74 ,"Tungsten     ","W ",183.84 ,110,0, 1.62],
[75 ,"Rhenium      ","Re",186.207,111,0, 1.51],
[76 ,"Osmium       ","Os",190.23 ,114,0, 1.44],
[77 ,"Iridium      ","Ir",192.217,115,0, 1.41],
[78 ,"Platinum     ","Pt",195.084,117,0, 1.36],
[79 ,"Gold         ","Au",196.967,118,0, 1.36],
[80 ,"Mercury      ","Hg",200.59 ,121,0, 1.32],
[81 ,"Thallium     ","Tl",204.383,123,0, 1.45],
[82 ,"Lead         ","Pb",207.2  ,125,0, 1.46],
[83 ,"Bismuth      ","Bi",208.98 ,126,0, 1.48],
[84 ,"Polonium     ","Po",210    ,126,1, 1.40],
[85 ,"Astatine     ","At",210    ,125,1, 1.50],
[86 ,"Radon        ","Rn",222    ,136,1, 1.50],
[87 ,"Francium     ","Fr",223    ,136,1, 2.60],
[88 ,"Radium       ","Ra",226    ,138,1, 2.21],
[89 ,"Actinium     ","Ac",227    ,138,1, 2.15],
[90 ,"Thorium      ","Th",232.038,142,1, 2.06],
[91 ,"Protactinium ","Pa",231.036,140,1, 2.00],
[92 ,"Uranium      ","U ",238.029,146,1, 1.96],
[93 ,"Neptunium    ","Np",237    ,144,1, 1.90],
[94 ,"Plutonium    ","Pu",244    ,150,1, 1.87],
[95 ,"Americium    ","Am",243    ,148,1, 1.80],
[96 ,"Curium       ","Cm",247    ,151,1, 1.69],
[97 ,"Berkelium    ","Bk",247    ,150,1, 2.00],
[98 ,"Californium  ","Cf",251    ,153,1, 2.00],
[99 ,"Einsteinium  ","Es",252    ,153,1, 2.00],
[100,"Fermium      ","Fm",257    ,157,1, 2.00],
[101,"Mendelevium  ","Md",258    ,157,1, 2.00],
[102,"Nobelium     ","No",259    ,157,1, 2.00],
[103,"Lawrencium   ","Lr",262    ,159,1, 2.00],
[104,"Rutherfordium","Rf",261    ,157,1, 2.00],
[105,"Dubnium      ","Db",262    ,157,1, 2.00],
[106,"Seaborgium   ","Sg",266    ,160,1, 2.00],
[107,"Bohrium      ","Bh",264    ,157,1, 2.00],
[108,"Hassium      ","Hs",267    ,159,1, 2.00],
[109,"Meitnerium   ","Mt",268    ,159,1, 2.00],
[110,"Darmstadtium ","Ds",271    ,161,1, 2.00],
[111,"Roentgenium  ","Rg",272    ,161,1, 2.00],
[112,"Copernicium  ","Cn",285    ,173,1, 2.00],
[113,"Nihonium     ","Nh",284    ,171,1, 2.00],
[114,"Flerovium    ","Fl",289    ,175,1, 2.00],
[115,"Moscovium    ","Mc",288    ,173,1, 2.00],
[116,"Livermorium  ","Lv",292    ,176,1, 2.00],
[117,"Tennessine   ","Ts",295    ,178,1, 2.00],
[118,"Oganesson    ","Og",294    ,176,1, 2.00],
]

# Create the elements object
elements = {}
for dat in data:
    d = {}
    for i in range(0, len(titles)):
        if isinstance(dat[i], str):    dat[i] = dat[i].strip() # Strip strings in dat
        if titles[i] == "radioactive": dat[i] = (dat[i] > 0)   # Replace radioactive with bool
        d[titles[i]] = dat[i]
    elements[dat[2]] = d

# Atom-substitution counts based on the ICSD 
# (useful for alchemical optimization)
#
# from 
#      https://tddft.org/bmg/files/data/pettifor/raw_data/substitution.dat
#
# based on the paper
#      https://doi.org/10.1088%2F1367-2630%2F18%2F9%2F093011
#
atom_substitutions = {}

# Parse from substitution.dat
with open(base_dir+"/substitution.dat") as f:
    i = -1
    for line in f:
        i += 1
        if i < 1: continue # First row/colum is padded with zeros
        ints = [int(w) for w in line.split()]

        # Get substitute atoms, sorted by frequency
        subs = {data[j-1][2] : c for j, c in enumerate(ints) if c > 0 and j > 0}
        subs = {k: v for k, v in sorted(subs.items(), key=lambda item: -item[1])}

        atom_substitutions[data[i-1][2]] = subs

def plot():
    import matplotlib.pyplot as plt
    import numpy as np
    
    an = []
    rs = []
    nc = []
    for e in elements:
        an.append(elements[e]["atomic number"])
        rs.append(elements[e]["covalent radius"])
        nc.append(elements[e]["neutrons"])

    p1 = plt.subplot(221)
    p2 = plt.subplot(222)
    p3 = plt.subplot(223)

    p1.plot(an,rs)
    p1.set_xlabel("Atomic number")
    p1.set_ylabel("Covalent radius (angstrom)")

    p2.plot(an,nc)
    p2.set_xlabel("Atomic number")
    p2.set_ylabel("Neutron count")

    mat = np.zeros((len(data), len(data)))
    for i in range(len(data)):
        for j in range(i):
            try: mat[i][j] = atom_substitutions[data[i][2]][data[j][2]]
            except: pass
    p3.imshow(mat)

    plt.show()
