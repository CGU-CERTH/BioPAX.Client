BioPAX.Client
=============

BioPAXClient is a Java application that can extract information from BioPAX files as well as create modified versions of template BioPAX files. The modified BioPAX files can contain presence or absence flags at their <i>Pathway</i> and/or <i>ProteinReference</i> elements through custom comments that follow a predefined format.


<h3>Build notes</h3>

- The project consists of a single file (<b>Main.java</b>) and requires an external library (<b>paxtools-4.1.1.jar</b>) in order to build it.
- The output of the building process is a JAR file (<b>BiopaxClient.jar</b> by default).
- Further instructions for distributing and running the  application are contained in the <b>/dist/README.txt</b> file.


<h3>Input parameters and Functionality</h3>

<table>
    <tr>
        <td align="center"><b>Input arguments</b></td>
        <td align="center"><b>Function</b></td>
    </tr>
    <tr>
        <td><b>-n</b> <i>[biopaxFile_input]</i></td>
        <td>Get pathway_organism name.</td>
    </tr>
    <tr>
        <td><b>-g</b> <i>[biopaxFile_input] ([genesList_output])</i></td>
        <td>Get genes list</td>
    </tr>
    <tr>
        <td><b>-i</b> <i>[biopaxFile_input] ([geneIds_output])</i></td>
        <td>Get genes custom ids</td>
    </tr>
    <tr>
        <td><b>-a</b> <i>[biopaxFile_input] [1/0]</i></td>
        <td>Add existence mark to the pathway</td>
    </tr>
    <tr>
        <td><b>-ap</b> <i>[biopaxFile_input] [protein_id] [1/0]</i></td>
        <td>Add existence mark to a protein of the pathway</td>
    </tr>
    <tr>
        <td><b>-an</b> <i>[biopaxFile_input] [pathNodeDir_path]</i></td>
        <td>Update pathway existence flag in a biopax file from node files information</td>
    </tr>
    <tr>
        <td><b>-agn</b> <i>[biopaxFile_input] [geneNodesDir_path]</i></td>
        <td>Update proteins existence flag in a biopax file from node files information</td>
    </tr>
    <tr>
        <td><b>-apgn</b> <i>[biopaxFile_input] [pathwayNodesDir_path] [genesNodeDir_path]</i></td>
        <td>Update pathway and proteins existence flag in a biopax file from node files information</td>
    </tr>
    <tr>
        <td><b>-c</b> <i>[biopaxFile_input]</i></td>
        <td>Check existence mark of the pathway</td>
    </tr>
    <tr>
        <td><b>-cp</b> <i>[biopaxFile_input] [protein_id]</i></td>
        <td>Check existence mark of a protein of the specified pathway</td>
    </tr>
    <tr>
        <td><b>-v</b> <i>[biopaxFile_input] [validation_output] [xml/html]</i></td>
        <td>Validate the biopax file</td>
    </tr>
    <tr>
        <td><b>-h</b></td>
        <td>Print help</td>
    </tr>
</table>

<br/>
Copyright (c) 2014 CERTH<br/>
Author: Dimitrios Vitsios<br/>
Last edit: 26 April 2014
