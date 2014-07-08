
void remove_spaces (std::string & fname)
{
	while (fname[0] == ' ') {
		fname = fname.substr(1);
	}
	while (fname[fname.size()-1] == ' ') {
		fname = fname.substr(0, fname.size() - 1);
	}
}


// -----------------------------------------------------------------------
//                               READ_SETS
// -----------------------------------------------------------------------
void read_sets (const std::string & file_name, std::map <std::string, std::vector <std::string> > & file_names, std::map <std::string, std::vector <std::string> > & bv_names)
{
	file_names.clear();
	bv_names.clear();
	
	std::ifstream infile;
	infile.open(file_name.c_str());
	if (!infile.good()) {
		std::cerr << "Cannot read file " << file_name << "\n";
		exit(1);
	}
	int nb_sets = 0;
	while (infile.good()) {
		std::string line;
		getline(infile, line);
		if (!line.empty()) {
			nb_sets++;
			std::stringstream current_tag;
			if (line.find(":") < line.size()) {
				current_tag << line.substr(0, line.find(":"));
				line = line.substr(line.find(":") + 1);
			} else {
				current_tag << "SET" << nb_sets;
			}
			std::vector <std::string> current_files;
			std::vector <std::string> current_bvs;
			while (!line.empty() && line.find(";") < line.size()) {
				std::string fname = line.substr(0, line.find(";"));
				remove_spaces(fname);
				std::string bv_name;
				if (fname.find(",") < fname.size()) {
					bv_name = fname.substr(fname.find(",") + 1);
					remove_spaces(bv_name);
					fname = fname.substr(0,fname.find(","));
					remove_spaces(fname);
				}
				current_files.push_back(fname);
				current_bvs.push_back(bv_name);
				line = line.substr(line.find(";") + 1);
			}
			std::string fname = line;
			remove_spaces(fname);
			std::string bv_name;
			if (fname.find(",") < fname.size()) {
				bv_name = fname.substr(fname.find(",") + 1);
				remove_spaces(bv_name);
				fname = fname.substr(0,fname.find(","));
				remove_spaces(fname);
			}
			current_files.push_back(fname);
			current_bvs.push_back(bv_name);
			file_names[current_tag.str()] = current_files;
			bv_names[current_tag.str()] = current_bvs;
		}
	}
	infile.close();
}
