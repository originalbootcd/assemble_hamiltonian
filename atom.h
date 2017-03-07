#ifndef ATOM_H
#define ATOM_H


class Atom
{
public:
	Atom() {
		neighbours = NULL;
		second_neighbours = NULL;
		x=0; y=0; z=0;
		dx=0; dy=0; dz=0;
		sostav=0;
	}

    ~Atom() {
		/*
		if (neighbours)
			delete neighbours;
		*/
        while (neighbours)
            Del(&neighbours);
		while (second_neighbours)
			Del(&second_neighbours);
    }

	double x, y, z;
	double dx, dy, dz;
    int sostav;
	List<Atom> *neighbours;
	List<Atom> *second_neighbours;
	int n;
};


class WaveFunction
{
public:
	WaveFunction() {
		atoms = NULL;
		atoms_tail = NULL;
	}

	~WaveFunction(){
		List<Atom> *head = atoms->next;
		while (head) {
			atoms->next = NULL;
			if (atoms)
				delete atoms;
			atoms = head;
			head = atoms->next;
		}
		if (atoms)
			delete atoms;
		/*
		if (atoms) {
			delete atoms;
			atoms = NULL;
		}*/
	}

	List<Atom> *atoms;
	List<Atom> *atoms_tail;
};




#endif
