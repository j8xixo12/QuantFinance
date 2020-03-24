#ifndef LATTICE_HPP_
#define LATTICE_HPP_

#include <vector>
#include <iostream>
#include <functional>

namespace LatticeMechanism {

    template <class Node, int LatticeType>
    class Lattice {
        private:
            std::vector<std::vector<Node>> tree;
            int type;

        public:
            Lattice() = delete;
            Lattice(const int& depth, const Node& val) {
                tree.resize(depth);
                int num = 1;
                for (auto& V_i : tree) {
                    V_i.reserve(num);
                    V_i.assign(num, val);
                    num++;
                }
            }
            Lattice(const Lattice<Node, LatticeType>& other) {
                tree.resize(other.depth);
                int num = 1;
                for (auto& V_i : tree) {
                    V_i.reserve(num);
                    V_i.assign(num, other.val);
                    num++;
                } 
            }
            virtual ~Lattice() {
                for (auto& V_i : tree) {
                    V_i.clear();
                }
                tree.clear();
            }

            const std::size_t Depth() const { return tree.capacity(); };
            std::vector<std::vector<Node>>& get_tree() { return tree; }

            Lattice<Node, LatticeType>& operator= (const Lattice<Node, LatticeType>& other) {
                tree.resize(other.depth);
                int num = 1;
                for (auto& V_i : tree) {
                    V_i.reserve(num);
                    V_i.assign(num, other.val);
                    num++;
                } 
            }

            std::vector<Node>& operator[](const int& index) { return tree[index]; }
    };

    template<class Node, int LatticeType>
    void print(Lattice<Node, LatticeType> &input) {
        for (const auto &Vector_i : input.get_tree()) {
            for (const auto &Vector_j : Vector_i) {
                std::cout << Vector_j << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    template<class Node, int LatticeType>
    void ModifyBaseValue(Lattice<Node, LatticeType> &lattice,
                            const std::function<Node (const Node& node)>& fun) {
        
        for (auto &Vector_i : lattice.get_tree()) {
            for (auto &Vector_j : Vector_i) {
                Vector_j = fun(Vector_j);
            }
        }
        return;                       
    }

    template<class Node, int LatticeType>
    void ForwardInduction(Lattice<Node, LatticeType> &lattice,
                            const std::function<std::tuple<Node, Node>(const Node& node)>& generator,
                            const Node& rootval) {
        lattice[0][0] = rootval;

        std::tuple<Node, Node> tuple;
        std::size_t v_size = 0;

        for (std::size_t i = 1; i < lattice.Depth(); ++i) {
            v_size = lattice[i - 1].capacity();
            for (std::size_t j = 0; j < v_size; ++j) {
                tuple = generator(lattice[i-1][j]);
                lattice[i][j] = std::get<0>(tuple);
                lattice[i][j + 1] = std::get<1>(tuple);
            }
        }
    }

    
    template<class Node, int LatticeType>
    Node BackwardInduction(Lattice<Node, LatticeType> &lattice,
                            const std::function<Node (const Node& upper, const Node& lower)>& generator,
                            const std::function<Node (const Node& node)>& endCondition) {
        auto it = lattice[lattice.Depth() - 1];
        for (auto i = it.begin(); i != it.end(); ++i) {
            (*i) = endCondition((*i));
        }

        for (std::size_t i = lattice.Depth() - 2; static_cast<int>(i) >= 0; --i) {
            for (std::size_t j = 0; j < lattice[i].capacity(); ++j) {
                lattice[i][j] = generator(lattice[i + 1][j + 1], lattice[i + 1][j]);
            }
        }

        return lattice[0][0];
    }
    
    template<class Node, int LatticeType>
    Node BackwardInduction(Lattice<Node, LatticeType> &lattice,
                            Lattice<Node, LatticeType> &l2,
                            const std::function<Node (const Node& upper, const Node& lower)>& generator,
                            const std::function<Node (const Node& node)>& endCondition,
                            const std::function<void (Node &node, const Node& val)>& constraintAdjuster) {        
        ModifyBaseValue(l2, endCondition);
        double tmp = 0.0;
        for (std::size_t i = l2.Depth() - 2; static_cast<int> (i) >= 0; --i) {
            for (std::size_t j = 0; j < l2[i].capacity(); ++j) {
                tmp = lattice[i][j];
                l2[i][j] = generator(l2[i + 1][j + 1], l2[i + 1][j]);
                constraintAdjuster(l2[i][j], tmp);
            }
        }
        return l2[0][0];
    }

    template<class Node, int LatticeType>
    Node BackwardInduction(Lattice<Node, LatticeType> &lattice,
                            Lattice<Node, LatticeType> &l2,
                            const std::function<Node (const Node& upper, const Node& lower)>& generator,
                            const std::function<Node (const Node& node)>& endCondition) {        
        ModifyBaseValue(l2, endCondition);
        double tmp = 0.0;
        for (std::size_t i = l2.Depth() - 2; static_cast<int> (i) >= 0; --i) {
            for (std::size_t j = 0; j < l2[i].capacity(); ++j) {
                tmp = lattice[i][j];
                l2[i][j] = generator(l2[i + 1][j + 1], l2[i + 1][j]);
            }
        }
        return l2[0][0];
    }
};

#endif // LATTICE_HPP_