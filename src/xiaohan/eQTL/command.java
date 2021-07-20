package xiaohan.eQTL;

import java.io.File;
import java.util.concurrent.Callable;

/**
 * @ author: yxh
 * @ created: 2021-07-20 : 10:43 PM
 */
class Command implements Callable<Command> {
    String command = null;
    File dir = null;

    public Command(String command, File dir) {
        this.command = command;
        this.dir = dir;
    }

    @Override
    public Command call() throws Exception {
        try {
            System.out.println(command);
            String[] cmdarry1 = {"/bin/bash", "-c", command};
            Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
            p1.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return this;
    }
}


