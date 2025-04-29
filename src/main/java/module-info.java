module solution.latticeboltzmannmodel.simulation.latticeboltzmannmodel {
    requires javafx.controls;
    requires javafx.fxml;


    opens solution.latticeboltzmannmodel.simulation to javafx.fxml;
    exports solution.latticeboltzmannmodel.simulation;
}