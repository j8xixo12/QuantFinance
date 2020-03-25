#ifndef LATTICE_HPP_
#define LATTICE_HPP_

#include <vector>
#include <iostream>
#include <functional>

namespace LatticeMechanism {

    template <class Node, int LatticeType>
    class Lattice {
        private:
            mutable std::vector<std::vector<Node>> tree;
            int type;

        public:
            Lattice() = delete;
            Lattice(const int& depth) {
                tree.resize(depth + 1);
                int num = 1;
                for(auto& V_i : tree) {
                    V_i.resize(num);
                    num++;
                }
            }
            Lattice(const int& depth, const Node& val) {
                tree.resize(depth + 1);
                int num = 1;
                for (auto& V_i : tree) {
                    V_i.reserve(num);
                    V_i.assign(num, val);
                    num++;
                }
            }
            Lattice(const Lattice<Node, LatticeType>& other) {
                tree.resize(other.Depth());
                int num = 1;
                for (auto& V_i : tree) {
                    V_i.resize(num);
                    num++;
                }
                std::copy(other.get_tree().begin(), other.get_tree().end(), tree.begin());
            }
            virtual ~Lattice() {
                for (auto& V_i : tree) {
                    V_i.clear();
                }
                tree.clear();
            }

            const std::size_t Depth() const { return tree.capacity(); };
            std::vector<std::vector<Node>>& get_tree() const { return tree; }

            Lattice<Node, LatticeType>& operator= (const Lattice<Node, LatticeType>& other) {
                tree.resize(other.Depth());
                int num = 1;
                for (auto& V_i : tree) {
                    V_i.resize(num);
                    num++;
                }
                std::copy(other.get_tree().begin(), other.get_tree().end(), tree.begin());
            }

            std::vector<Node>& operator[](const int& index) const { return tree[index]; }
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
    void ModifyFinalValue(Lattice<Node, LatticeType> &lattice,
                            const std::function<Node (const Node& node)>& fun) {
        for (auto &i : lattice[lattice.Depth() - 1]) {
            i = fun(i);
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
        ModifyFinalValue(lattice, endCondition);

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
        ModifyFinalValue(l2, endCondition);
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
    Lattice<std::tuple<Node, Node>, LatticeType> merge (const Lattice<Node, LatticeType>& latticeA,
                                                        const Lattice<Node, LatticeType>& latticeB) {
        Lattice<std::tuple<Node, Node>, LatticeType> result(latticeA.Depth());
        
        for(auto i = 0; i < latticeA.Depth(); ++i) {
            for(auto j = 0; j < (latticeA[i]).capacity(); ++j) {
                result[i][j] = std::make_tuple(latticeA[i][j], latticeB[i][j]);
            }
        }
        return result;
    }
};

#endif // LATTICE_HPP_