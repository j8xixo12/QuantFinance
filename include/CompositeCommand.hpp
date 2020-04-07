#ifndef COMPOSITECOMMAND_HPP_
#define COMPOSITECOMMAND_HPP_

#include <list>
#include "OptionCommand.hpp"

using CompositeCommandData = std::vector<std::shared_ptr<OptionCommand>>;

class CompositeCommand : public OptionCommand {
    private:
        CompositeCommandData data;
        std::vector<double> v_data;
    public:
        CompositeCommand() : data(CompositeCommandData()) {}

        void add(OptionCommand* cmd) {
            data.push_back(std::shared_ptr<OptionCommand>(cmd)); 
        }
        
        std::size_t size() const {
            return data.size(); 
        }

        double execute(double S) override { // Iterate over each element of the composite and call 
            // its execute()
            v_data.resize(data.size());
            for (auto i = 0; i < data.size(); ++i) {
                v_data[i] = data[i]->execute(S); // nested pointer 
            }
        }
};

#endif // COMPOSITECOMMAND_HPP_