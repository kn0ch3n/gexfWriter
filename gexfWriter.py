import xml.etree.ElementTree as ET
import xml.dom.minidom
import os
import csv
from collections import defaultdict, OrderedDict

### class for creating a gexf file
class Gexf(object):
	def __init__(self, sbmlversion = "{http://www.sbml.org/sbml/level2/version3}", modelfile = None, metabolitefiles = [], genefiles = [], genelist = [],
				compartment = None, separate_genes = False,
				attribs = {"xmlns": "http://www.gexf.net/1.2draft", "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
							"xsi:schemaLocation": "http://www.gexf.net/1.2draft http://www.gexf.net/1.2draft/gexf.xsd", "version": "1.2"}):
		"""
		Creates a gexf file, including xml & namespace declarations, and an empty graph tag.
		Also creates these variables:
			- modelfile
			- modelroot
			- speciesroot
			- reactionsroot
			- compartment (e.g. "C_c" for cytosol)
		"""
		self.separate_genes = separate_genes # if a gene acts in multiple reactions, should there be separate instances or only one for all? TODO: make this better/rename/something like that
		self.sbmlversion = sbmlversion
		self.modelfile = modelfile
		self.metabolitefiles = metabolitefiles
		self.genefiles = genefiles
		self._genelist = genelist
		self._genefilters = {} # TODO: best practice?
		self._genefilteroperator = 'or'
		self.compartment = compartment # TODO: maybe remove this here, not sure what's the best way of handling the compartments...
		model = ET.parse(self.modelfile)
		self.modelroot = model.getroot()
		self.speciesroot = self.modelroot.find("./{0}model/{0}listOfSpecies".format(self.sbmlversion))
		self.reactionsroot = self.modelroot.find("./{0}model/{0}listOfReactions".format(self.sbmlversion))
		self.root = ET.Element("gexf")
		for k, v in attribs.items():
			self.root.attrib[k] = v
		self.graph = self.build_node(self.root, "graph")
	
	# The xml file of the model
	# TODO: documentation? how?
	@property
	def modelfile(self):
		return self._modelfile
	@modelfile.setter
	def modelfile(self, value):
		self._modelfile = value
	
	# A list of 2-tuples: [i][0] is the prefix that is given to the columns in the value file. [i][1] is the name of a csv file with values. TODO: update documentation
	@property
	def metabolitefiles(self):
		return self._metabolitefiles
	@metabolitefiles.setter
	def metabolitefiles(self, value):
		self._metabolitefiles = value
	
	@property
	def genefiles(self):
		return self._genefiles
	@genefiles.setter
	def genefiles(self, value):
		self._genefiles = value
	
	# allows filtering genes by their name
	# only genes in the list will be included
	@property
	def genelist(self):
		return self._genelist
	@genelist.setter
	def genelist(self, value):
		self._genelist = value
	
	# allows filtering genes by values of their attributes
	# is a dict: key is the "for" value that will be filtered
	# value is a function that takes one input argument (the value the key refers to) and returns True if it passes of False if the value is filtered
	@property
	def genefilters(self):
		return self._genefilters
	@genefilters.setter
	def genefilters(self, value):
		self._genefilters = value
	
	# you can add multiple gene filters, then a gene must pass only one of them by default ( = "or"), e.g. 
	# you can set this to "and", then a gene must pass every filter in genefilters to pass
	@property
	def genefilteroperator(self):
		return self._genefilteroperator
	@genefilteroperator.setter
	def genefilteroperator(self, value):
		value = str(value).lower()
		if not value in ["or", "and"]:
			raise IOError("genefilteroperator must be either 'or' or 'and'. Got: " + str(value))
		else:
			self._genefilteroperator = value
	
	def node_attributes(self):
		"""
		Returns a list of OrderedDicts, where each OrderedDict represents a node_attribute, e.g.:
			<attribute id="0" title="type" type="string"/>
		would be:
			[OrderedDict([('id', 0), ('title', 'type'), ('type', 'string')])]
		
		The above example is always the first attribute with the id "0".
		The other attributes are read from the first line of metabolitefiles and genefiles,
		thus corresponding to the column names of the node (aka gene or metabolite) values.
		"""
		# TODO: make this extendable? like using types other than float, allow only certain values, etc...
		result = []
		# add the node type as first attribute
		result.append(OrderedDict({"id": "0", "title": "type", "type": "string"}))
		# TODO: allow only "metabolite", "gene" and "reaction" for type
		nid = 1 # node id
		# metabolite attributes & gene attributes
		for files in [self.metabolitefiles, self.genefiles]:
			for d in files:
				if not "filename" in d:
					raise IOError("No filename given...")
				# read the first line of every metabolitefile
				with open(d["filename"], 'r') as f:
					reader = csv.reader(f, delimiter='\t', quotechar='"')
					colnames = next(reader)[1:] # skip the first column, for it is empty
				for name in colnames:
					result.append(OrderedDict({"id": str(nid), "title": d["prefix"] + name if "prefix" in d else name, "type": "float"}))
					nid += 1
		return result
	
	def attvalues(self):
		"""
		Sets self.attvalues to a dict where the key is a node id, and the value is a dict.
		e.g.:
			<node id="M_m0417" label="[protein]-L-lysine (C_n)">
			<attvalues>
				<attvalue for="0" value="4.40528634361"/>
				<attvalue for="1" value="2.40528634361"/>
			</attvalues>
		</node>
		would be:
			{"M_m0417": {"0": 4.40528634361, "1": 2.40528634361}}
		"""
		result = defaultdict(dict)
		startid = 1
		# metabolite values & gene values
		for files in [self.metabolitefiles, self.genefiles]:
			for d in files:
				if not "filename" in d:
					raise IOError("No filename given...")
				lines = []
				with open(d["filename"], 'r') as f:
					reader = csv.reader(f, delimiter='\t', quotechar='"')
					columnlength = len(next(reader)) - 1 # remove the header, and keep the number of values in each line
					for l in reader:
						lines.append(l)
				for l in lines:
					id = startid
					label = l[0] # note that metabolites have the compartment here for indexing, genes don't
					values = l[1:]
					for v in values:
						result[label].update({str(id): str(v)})
						id += 1
				startid += columnlength
		return result
	
	def set_attvalues(self):
		self.attvalues = self.attvalues()
	
	def build(self):
		"""
		This actually builds the xml structure by modifying self.root.
		"""
		self.build_node_attributes()
		self.set_attvalues()
		self.build_nodes_and_edges()
	
	def build_node_attributes(self):
		"""
		Adds node attributes to the graph, calling node_attributes() to get them. Is called by build.
		e.g.:
		<attributes class="node">
			<attribute id="0" title="Female_Non-directional" type="float"/>
		</attributes>
		"""
		attributes = self.build_node(self.graph, "attributes")
		attributes.attrib["class"] = "node"
		for od in self.node_attributes():
			attribute = self.build_node(attributes, "attribute")
			for k, v in od.items():
				attribute.attrib[str(k)] = str(v)
				
	def build_gene(self, parent, gid, label, geneattribs):
		node = self.build_node(self.nodes, "node", {"id": gid, "label": label}) # add gene node
		attvalues = self.build_node(node, "attvalues") # add attvalues to gene node
		# add gene attributes
		self.build_node(attvalues, "attvalue", {"for": "0", "value": "gene"})
		for k, v in geneattribs.items():
			self.build_node(attvalues, "attvalue", {"for": k, "value": v})
	
	def build_node(self, parent, name, attribs={}):
		"""
		Basically a wrapper for ElementTree.SubElement().
		Adds a node called name as a child to parent, with the attributes from attribs.
		Returns the node.
		"""
		node = ET.SubElement(parent, name)
		for k, v in attribs.items():
			node.attrib[str(k)] = str(v)
		return node
	
	def build_metabolites(self):
		"""
		Adds metabolite nodes to the DOM. Is called by build_nodes_and_edges().
		"""
		for m in self.metabolites:
			label = m.attrib["name"] # e.g. (10Z)-heptadecenoic acid
			node = self.build_node(self.nodes, "node", {"id": m.attrib["id"], "label": label})
			# add node attributes
			attvalues = self.build_node(node, "attvalues")
			# first the node type...
			self.build_node(attvalues, "attvalue", {"for": "0", "value": "metabolite"})
			# ... then the rest...
			for k, v in self.attvalues[label + " (" + compartment + ")"].items():
				self.build_node(attvalues, "attvalue", {"for": k, "value": v})
	
	def build_nodes_and_edges(self):
		"""
		Adds metabolite nodes to the graph. 
		Also saves a list of metabolites in this compartment and all genes.
		e.g.:
		<nodes>
		<node id="M_m0010" label="(10Z)-heptadecenoyl-CoA (C_m)">
			<attvalues>
				<attvalue for="0" value="2.63157894737"/>
			</attvalues>
		</node>
		</nodes>
		"""
		# create empty nodes and edges
		self.nodes = self.build_node(self.graph, "nodes")
		self.edges = self.build_node(self.graph, "edges")
		# start adding stuff
		self.species = self.speciesroot.findall("./*") # TODO: one could filter here with xpath, but it makes my head hurt...
		# all metabolite ids start with "M_m", and we only want those in our compartment
		self.metabolites = [s for s in self.species if s.attrib['compartment'] == self.compartment and s.attrib['id'].startswith("M_m")]
		# all gene ids start with "E_", and they are all in C_c, so the compartment does not matter
		self.genes = [s for s in self.species if s.attrib['id'].startswith("E_")]
		# add metabolite nodes
		self.build_metabolites()
		# add gene nodes and edges
		self.build_genes_and_reactions()
	
	def build_genes_and_reactions(self):
		"""
		Adds gene nodes to the DOM. Also adds all the edges (aka reactions).
		Would also add reaction nodes if those were still implemented.
		"""
		metabolite_ids = set([m.attrib["id"] for m in self.metabolites])
		reactions = self.reactionsroot.findall("./*")
		genes = {} # keep track of the genes we have added already: key is label, value is the node
		eid = 0 # edge id
		if self.separate_genes: # we can't use the gene identifiers from the model, so use gid
			gid = 0 # gene id
		for r in reactions:
			# get reactants, products and modifier genes for that reaction
			reactants = r.findall("./{0}listOfReactants/{0}speciesReference".format(self.sbmlversion))
			products = r.findall("./{0}listOfProducts/{0}speciesReference".format(self.sbmlversion))
			reactiongenes = r.findall("./{0}listOfModifiers/{0}modifierSpeciesReference".format(self.sbmlversion))
			# only add the reaction & associated genes if all reactants and products are in this compartment
			if (set([m.attrib["species"] for m in reactants]).issubset(metabolite_ids)):
				#and set([m.attrib["species"] for m in products]).issubset(metabolite_ids)): # if this line is commented out, don't filter by products
				# add genes for that reaction
				for gene in reactiongenes:
					if self.genefilteroperator == "or" and len(self.genefilters) > 0:
						gene_passes = False
					elif self.genefilteroperator == "and" or len(self.genefilters) == 0:
						gene_passes = True # set it to true True, because with "and" a gene is skipped if it does not pass any of the filters
					# get a label for the gene
					label = gene.attrib["species"] # E_xxx
					for s in self.species: # TODO: not sure if this is necessary, and it's probably a bit slow
						if s.attrib["id"] == gene.attrib["species"]:
							label = s.attrib["name"] # ENSG_xxx
					# skip the gene if it is not in the genelist
					if self.genelist != [] and label not in self.genelist:
						continue
					# skip the gene if it does not pass the filters
					for k, v in self.genefilters.items():
						if not self.attvalues[label] or not v(self.attvalues[label][k]): # gene does not pass the filter
							if self.genefilteroperator == "and": # only if the gene must pass all filters (aka boolean AND)...
								continue # ... don't add this gene
						else: # gene passes the filter
							if self.genefilteroperator == "or":
								gene_passes = True
					else: # if we didn't continue, add gene as a node and reactions as edges
						if not gene_passes:
							continue # don't add this gene
						if not self.separate_genes: # make only one node for a gene, even if it acts in multiple reactions
							gid = gene.attrib["species"]
							if gid in genes: # use the gene node if it already exists
								node = genes[gid] 
							else: # otherwise create a new gene node
								node = self.build_gene(self.nodes, gid, label, self.attvalues[label])
						else: # we want separate genes, so add a new gene node anyway
							node = self.build_gene(self.nodes, gid, label, self.attvalues[label])
						# add edge from each reactant to gene
						for reactant in reactants:
							self.build_node(self.edges, "edge", {"id": eid, "source": reactant.attrib["species"], "target": gid})
							eid += 1
						# add edge from gene to each product
						for product in products:
							self.build_node(self.edges, "edge", {"id": eid, "source": gid, "target": product.attrib["species"]})
							eid += 1
						if self.separate_genes:
							gid += 1
						genes[gid] = node
		
	def pretty(self):
		"""
		Returns the current xml content with indentations and linebreaks.
		"""
		return xml.dom.minidom.parseString(ET.tostring(self.root)).toprettyxml()
		
	def write(self, filename = None):
		"""
		Writes the xml to a file. If no filename is given: <modelfile without extension>.gexf.
		"""
		if filename == None:
			filename = self.modelfile.split(".")[0] + ".gexf"
		with open(filename, "w") as f:
			f.write(self.pretty())
			print("File written to " + os.path.abspath(filename))
			
	def __str__(self):
		"""
		Returns the current xml content, without indentations or linebreaks. Use pretty() if you want those.
		"""
		result = ET.tostring(self.root)
		return str(result)

# genes for "Chondroitin sulfate biosynthesis"
genes = ["ENSG00000114646", "ENSG00000132692", "ENSG00000182492", "ENSG00000130287", "ENSG00000038427", "ENSG00000173546", "ENSG00000011465", "ENSG00000182022", "ENSG00000147408", "ENSG00000169826", "ENSG00000123989", "ENSG00000033100", "ENSG00000131873", "ENSG00000198108", "ENSG00000171310", "ENSG00000136213", "ENSG00000180767", "ENSG00000154080", "ENSG00000122863", "ENSG00000147119"]
compartments = ["C_s", "C_p", "C_m", "C_c", "C_l", "C_r", "C_g", "C_n"]
gexfs = {}

for compartment in compartments:
	print("Starting " + compartment + "...")
	gexfs[compartment] = Gexf(modelfile="input/iAdipocytes1809.xml", compartment = compartment, separate_genes=False)
	gexfs[compartment].metabolitefiles = [{"prefix": "Female_", "filename": "input/metabolite_reporter_consrank_female.csv"}, {"prefix": "Male_", "filename": "input/metabolite_reporter_consrank_male.csv"}] # csv files with statistics for metabolites
	gexfs[compartment].genefiles = [{"prefix": "Female_", "filename": "input/pfc_ENSG_female.csv"}, {"prefix": "Male_", "filename": "input/pfc_ENSG_male.csv"}]
	#gexfs[compartment].genefilters = {"21": lambda x: float(x) > 20, "23": lambda x: float(x) > 20} # only female gene p-value < 0.05 ( = inverted p-value > 20) OR male gene p-value < 0.05
	gexfs[compartment].genelist = genes
	gexfs[compartment].build() # builds the xml from the given files... might put this into write()
	gexfs[compartment].write(filename = "output/" + compartment + "_chondroitin_sulfate_biosynthesis.gexf")
