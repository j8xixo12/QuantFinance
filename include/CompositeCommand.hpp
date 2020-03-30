#ifndef COMPOSITECOMMAND_HPP_
#define COMPOSITECOMMAND_HPP_

#include <list>
#include "OptionCommand.hpp"

using CompositeCommandData = std::list<std::shared_ptr<OptionCommand>>;

class CompositeCommand : public OptionCommand {
    private:
        CompositeCommandData data;
    public:
        CompositeCommand() : data(CompositeCommandData()) {}

        void add(OptionCommand* cmd) {
            data.push_back(std::shared_ptr<OptionCommand>(cmd)); 
        }
        
        std::size_t size() const {
            return data.size(); 
        }

        void execute(double S) override { // Iterate over each element of the composite and call 
            // its execute()
            for (auto it = data.begin(); it != data.end(); ++it) {
                (*it)->execute(S); // nested pointer 
            }
        }
};

#endif // COMPOSITECOMMAND_HPP_