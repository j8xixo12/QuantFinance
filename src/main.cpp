#include "OptionCommand.hpp"
#include "CompositeCommand.hpp"



void TestCP() { // Haug 2007 page 3, C
    double K = 65.0; 
    double T = 0.25; 
    double r = 0.08; 
    double b = r; 
    double sig = 0.3;
    double S = 60.0;

    CallPrice cp(K, T, r, b, sig); 
    cp.execute(S);
}


void TestCompositeCommand() { // Create a composite option command
    double K = 65.0; 
    double T = 0.25; 
    double r = 0.08; 
    double b = r; 
    double sig = 0.3;
    double S = 60.0;
    // Create composite command 
    std::cout << "Command on options\n"; 
    CompositeCommand optionCommand;
    // Add prices and their deltas 
    optionCommand.add(new CallPrice (K, T, r, b, sig)); 
    // optionCommand.add(new PutPrice(K, T, r, b, sig));
    optionCommand.add(new CallDelta(K, T, r, b, sig)); 
    // optionCommand.add(new PutDelta(K, T, r, b, sig));
    std::cout << "Size of commandlist: " << optionCommand.size();
    optionCommand.execute(S);
}

int main(int argc, char* argv[]) {
    TestCP();
    TestCompositeCommand();
    return 0;
}
