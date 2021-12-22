//
// This file holds functions used to get relative publish folder for igv.js.
// It is not a real solution. It will not handle the user defined configs.
//

class RelativePublishFolder {
    //
    // Check publishDir from configFiles
    //
    public static String getPublishedFolder(workflow, moduleName){
        def p = "${moduleName.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}"
        def parser = new ConfigSlurper()
        ConfigObject properties = parser.parse(new File("${workflow.projectDir}/conf/modules.config").text)
        def prop = properties.process[moduleName]
        if(prop){
            p = prop.publishDir?prop.publishDir.path?prop.publishDir.path().replace("[:]/",""):p:p
        }
        return(p+'/')
    }
}
