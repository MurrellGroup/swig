declare module "@biowasm/aioli" {
  interface MountedFile {
    name: string;
    data: string | File | Blob | Uint8Array;
  }

  interface AioliRuntime {
    mount(file: MountedFile): Promise<string>;
    exec(command: string): Promise<string>;
  }

  export default class Aioli {
    constructor(tools: string[], config?: Record<string, unknown>);
    then<TResult1 = AioliRuntime, TResult2 = never>(
      onfulfilled?: ((value: AioliRuntime) => TResult1 | PromiseLike<TResult1>) | null,
      onrejected?: ((reason: unknown) => TResult2 | PromiseLike<TResult2>) | null,
    ): Promise<TResult1 | TResult2>;
  }
}
